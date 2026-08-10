#!/usr/bin/env python3
"""校验 compute-molecular-features 输出契约和科学失败关闭规则。"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import sys
from pathlib import Path
from typing import Any


SCHEMA_VERSION = "1.0.0"
WORKFLOW = "molecular-feature-computation"
DESCRIPTOR_SET_ID = "rdkit-2d-core-v1"
CALCULATION_VIEWS = {"standardized", "parent"}
CALCULATION_STATUSES = {"completed", "partial", "not_run", "error"}
DISPOSITIONS = {"ready_for_downstream", "review_required", "rejected"}
TEMPORAL_KEYS = {"generated_at_utc", "retrieved_at_utc", "requested_at_utc"}
FINGERPRINT_NAMES = {"morgan", "rdkit_topological", "maccs"}
DESCRIPTOR_NAMES = {
    "MolecularFormula",
    "MolecularWeight",
    "ExactMolWt",
    "HeavyAtomCount",
    "NumHDonors",
    "NumHAcceptors",
    "NumRotatableBonds",
    "RingCount",
    "NumAromaticRings",
    "FractionCSP3",
    "TPSA",
    "MolLogP",
    "FormalCharge",
    "NumHeteroatoms",
}
REQUIRED_TOP_LEVEL = {
    "schema_version",
    "workflow",
    "generated_at_utc",
    "tool_versions",
    "options",
    "descriptor_set",
    "fingerprint_profiles",
    "input_summary",
    "records",
    "dataset_profile",
    "errors",
    "warnings",
    "notices",
    "human_review_required",
    "result_fingerprint",
}
REQUIRED_RECORD_FIELDS = {
    "id",
    "record_index",
    "original_structure",
    "standardized_structure",
    "parent_structure",
    "source_structure",
    "calculation_view",
    "calculation_status",
    "descriptors",
    "fingerprints",
    "missing_features",
    "qc_findings",
    "upstream_parse_status",
    "upstream_standardization_status",
    "upstream_disposition",
    "upstream_human_review_required",
    "upstream_workflow",
    "upstream_fingerprint",
    "upstream_tool_versions",
    "upstream_profile",
    "input_record_fingerprint",
    "disposition",
    "human_review_required",
}
SECRET_RE = re.compile(
    r"ark-[A-Za-z0-9_-]{12,}|"
    r"Bearer\s+[A-Za-z0-9._~+/=-]{12,}|"
    r"(?:Authorization|Cookie|Token|Api[_ -]?Key)\s*[:=]\s*\S{12,}",
    re.IGNORECASE,
)
FORBIDDEN_CLAIMS = {
    "药效已确认",
    "活性已确认",
    "毒性已确认",
    "安全性已确认",
    "结构已确证",
    "适合直接建模",
    "suitable for modeling",
    "same biological function",
    "proven active",
    "proven safe",
    "safe to synthesize",
}


def canonical_json(value: Any) -> str:
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"))


def sha256_text(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def sha256_json(value: Any) -> str:
    return sha256_text(canonical_json(value))


def _without_temporal_fields(value: Any) -> Any:
    if isinstance(value, dict):
        return {
            key: _without_temporal_fields(item)
            for key, item in value.items()
            if key not in TEMPORAL_KEYS and key != "result_fingerprint"
        }
    if isinstance(value, list):
        return [_without_temporal_fields(item) for item in value]
    return value


def output_fingerprint(document: dict[str, Any]) -> str:
    return sha256_json(_without_temporal_fields(document))


def require_keys(
    value: dict[str, Any],
    required: set[str],
    path: str,
    errors: list[str],
) -> None:
    missing = sorted(required - set(value))
    if missing:
        errors.append(f"{path} missing fields: {missing!r}")


def find_non_finite(value: Any, path: str = "$") -> list[str]:
    results = []
    if isinstance(value, float) and not math.isfinite(value):
        results.append(path)
    elif isinstance(value, dict):
        for key, item in value.items():
            results.extend(find_non_finite(item, f"{path}.{key}"))
    elif isinstance(value, list):
        for index, item in enumerate(value):
            results.extend(find_non_finite(item, f"{path}[{index}]"))
    return results


def validate_finding(item: Any, path: str, errors: list[str]) -> None:
    if not isinstance(item, dict):
        errors.append(f"{path} must be an object")
        return
    for field in ("code", "severity", "message", "source"):
        if not item.get(field):
            errors.append(f"{path}.{field} is required")
    if item.get("severity") not in {"error", "warning", "review", "notice"}:
        errors.append(f"{path}.severity is invalid")


def validate_profile(
    name: str, profile: Any, path: str, errors: list[str]
) -> None:
    if not isinstance(profile, dict):
        errors.append(f"{path} must be an object")
        return
    require_keys(
        profile,
        {
            "profile_id",
            "algorithm",
            "method_family",
            "representation",
            "parameters",
            "known_limitations",
            "profile_fingerprint",
        },
        path,
        errors,
    )
    if profile.get("representation") != "bit_vector_on_bits":
        errors.append(f"{path}.representation must be bit_vector_on_bits")
    expected_profile_fingerprint = sha256_json(
        {
            key: value
            for key, value in profile.items()
            if key != "profile_fingerprint"
        }
    )
    if profile.get("profile_fingerprint") != expected_profile_fingerprint:
        errors.append(f"{path}.profile_fingerprint mismatch")
    parameters = profile.get("parameters")
    if not isinstance(parameters, dict):
        errors.append(f"{path}.parameters must be an object")
        return
    if name == "morgan":
        for field in (
            "radius",
            "fpSize",
            "includeChirality",
            "useBondTypes",
            "countSimulation",
            "includeRedundantEnvironments",
            "bitsPerFeature",
        ):
            if field not in parameters:
                errors.append(f"{path}.parameters.{field} is required")
    elif name == "rdkit_topological":
        for field in (
            "minPath",
            "maxPath",
            "useHs",
            "branchedPaths",
            "useBondOrder",
            "countSimulation",
            "fpSize",
            "numBitsPerFeature",
        ):
            if field not in parameters:
                errors.append(f"{path}.parameters.{field} is required")
    elif name == "maccs":
        if parameters.get("fpSize") != 167:
            errors.append(f"{path}.parameters.fpSize must be 167")
        if parameters.get("bit0Unused") is not True:
            errors.append(f"{path}.parameters.bit0Unused must be true")


def validate_fingerprint(
    fingerprint: Any,
    profile: dict[str, Any],
    path: str,
    errors: list[str],
) -> None:
    if not isinstance(fingerprint, dict):
        errors.append(f"{path} must be an object")
        return
    require_keys(
        fingerprint,
        {
            "profile_id",
            "representation",
            "size",
            "on_bits",
            "bit_count",
            "density",
            "bitvector_sha256",
            "hash_encoding",
        },
        path,
        errors,
    )
    if fingerprint.get("profile_id") != profile.get("profile_id"):
        errors.append(f"{path}.profile_id does not match fingerprint profile")
    if fingerprint.get("representation") != "bit_vector_on_bits":
        errors.append(f"{path}.representation is invalid")
    size = fingerprint.get("size")
    on_bits = fingerprint.get("on_bits")
    if not isinstance(size, int) or size <= 0:
        errors.append(f"{path}.size must be a positive integer")
        return
    if not isinstance(on_bits, list) or not all(
        isinstance(item, int) for item in on_bits
    ):
        errors.append(f"{path}.on_bits must be an integer list")
        return
    if on_bits != sorted(set(on_bits)):
        errors.append(f"{path}.on_bits must be sorted and unique")
    if any(item < 0 or item >= size for item in on_bits):
        errors.append(f"{path}.on_bits contains an out-of-range bit")
    if fingerprint.get("bit_count") != len(on_bits):
        errors.append(f"{path}.bit_count does not match on_bits")
    expected_density = len(on_bits) / size
    density = fingerprint.get("density")
    if not isinstance(density, (int, float)) or not math.isclose(
        float(density), expected_density, rel_tol=0.0, abs_tol=1e-15
    ):
        errors.append(f"{path}.density does not match bit_count/size")
    bit_set = set(on_bits)
    ascii_bits = "".join(
        "1" if index in bit_set else "0" for index in range(size)
    )
    if fingerprint.get("bitvector_sha256") != sha256_text(ascii_bits):
        errors.append(f"{path}.bitvector_sha256 mismatch")
    if fingerprint.get("hash_encoding") != (
        "ascii_bitstring_index_0_to_n_minus_1"
    ):
        errors.append(f"{path}.hash_encoding is invalid")


def validate_record(
    record: Any,
    index: int,
    descriptor_names: set[str],
    profiles: dict[str, Any],
    calculation_view: str,
    errors: list[str],
    warnings: list[str],
) -> None:
    path = f"records[{index}]"
    if not isinstance(record, dict):
        errors.append(f"{path} must be an object")
        return
    require_keys(record, REQUIRED_RECORD_FIELDS, path, errors)
    if not REQUIRED_RECORD_FIELDS <= set(record):
        return
    if not isinstance(record["id"], str) or not record["id"]:
        errors.append(f"{path}.id must be a non-empty string")
    if record["record_index"] != index:
        errors.append(f"{path}.record_index must preserve input order")
    if not isinstance(record["original_structure"], str):
        errors.append(f"{path}.original_structure must be a string")
    if record["calculation_view"] != calculation_view:
        errors.append(f"{path}.calculation_view does not match options")
    if record["calculation_status"] not in CALCULATION_STATUSES:
        errors.append(f"{path}.calculation_status is invalid")
    if record["disposition"] not in DISPOSITIONS:
        errors.append(f"{path}.disposition is invalid")
    expected_source = (
        record["standardized_structure"]
        if calculation_view == "standardized"
        else record["parent_structure"]
    )
    if record["source_structure"] != expected_source:
        errors.append(f"{path}.source_structure mixes calculation views")
    for field in (
        "descriptors",
        "fingerprints",
    ):
        if not isinstance(record[field], dict):
            errors.append(f"{path}.{field} must be an object")
    for field in (
        "missing_features",
        "qc_findings",
        "upstream_human_review_required",
        "human_review_required",
    ):
        if not isinstance(record[field], list):
            errors.append(f"{path}.{field} must be a list")
    for finding_index, item in enumerate(record["qc_findings"]):
        validate_finding(
            item, f"{path}.qc_findings[{finding_index}]", errors
        )

    upstream_blocks = (
        record["upstream_disposition"] == "rejected"
        or record["upstream_parse_status"] == "error"
        or record["upstream_standardization_status"] in {"error", "not_run"}
    )
    if upstream_blocks:
        if record["calculation_status"] != "not_run":
            errors.append(f"{path} rejected upstream record must not run")
        if record["descriptors"] or record["fingerprints"]:
            errors.append(f"{path} rejected upstream record emitted features")
        if record["disposition"] != "rejected":
            errors.append(f"{path} rejected upstream record must stay rejected")

    upstream_review = {
        (
            str(item["code"])
            if isinstance(item, dict) and item.get("code")
            else sha256_json(item)
            if isinstance(item, dict)
            else str(item)
        )
        for item in record["upstream_human_review_required"]
    }
    if record["upstream_disposition"] == "review_required":
        if record["disposition"] == "ready_for_downstream":
            errors.append(f"{path} lost upstream review_required disposition")
    if upstream_review and not upstream_review <= set(
        record["human_review_required"]
    ):
        errors.append(f"{path} lost upstream human review reasons")

    status = record["calculation_status"]
    if status in {"completed", "partial"}:
        if not isinstance(record["source_structure"], str) or not record[
            "source_structure"
        ]:
            errors.append(f"{path} calculated without a source structure")
        if not isinstance(record["descriptors"], dict):
            errors.append(f"{path}.descriptors must be an object")
        elif set(record["descriptors"]) != descriptor_names:
            errors.append(f"{path}.descriptors does not match descriptor_set")
        if not isinstance(record["fingerprints"], dict):
            errors.append(f"{path}.fingerprints must be an object")
        elif set(record["fingerprints"]) != FINGERPRINT_NAMES:
            errors.append(
                f"{path}.fingerprints does not match fingerprint_profiles"
            )
        if isinstance(record["fingerprints"], dict):
            for name, profile in profiles.items():
                fingerprint = record["fingerprints"].get(name)
                if fingerprint is not None:
                    validate_fingerprint(
                        fingerprint,
                        profile,
                        f"{path}.fingerprints.{name}",
                        errors,
                    )
    if status == "completed" and record["missing_features"]:
        errors.append(f"{path} completed record cannot have missing_features")
    if status == "partial" and not record["missing_features"]:
        errors.append(f"{path} partial record must list missing_features")
    if status in {"not_run", "error"}:
        if record["descriptors"] or record["fingerprints"]:
            errors.append(f"{path} non-calculated record emitted features")
    if status == "error" and record["disposition"] != "rejected":
        errors.append(f"{path} calculation error must be rejected")
    if record["disposition"] == "ready_for_downstream":
        if status != "completed":
            errors.append(f"{path} ready record must be completed")
        if record["human_review_required"]:
            errors.append(f"{path} ready record cannot require human review")
    if record["disposition"] != "ready_for_downstream":
        warnings.append(f"{path} disposition is {record['disposition']}")


def validate_dataset_profile(
    profile: Any,
    records: list[dict[str, Any]],
    descriptor_names: set[str],
    errors: list[str],
) -> None:
    if not isinstance(profile, dict):
        errors.append("dataset_profile must be an object")
        return
    for field in (
        "total_records",
        "calculation_status_counts",
        "disposition_counts",
        "descriptor_statistics",
        "constant_features",
        "near_constant_features",
        "duplicate_structures",
        "fingerprint_density_statistics",
        "human_review_count",
        "statistical_qc_parameters",
        "interpretation",
    ):
        if field not in profile:
            errors.append(f"dataset_profile.{field} is required")
    if profile.get("total_records") != len(records):
        errors.append("dataset_profile.total_records does not match records")
    expected_statuses = {
        status: sum(
            record.get("calculation_status") == status for record in records
        )
        for status in CALCULATION_STATUSES
    }
    if profile.get("calculation_status_counts") != expected_statuses:
        errors.append(
            "dataset_profile.calculation_status_counts does not match records"
        )
    expected_dispositions = {
        status: sum(record.get("disposition") == status for record in records)
        for status in DISPOSITIONS
    }
    if profile.get("disposition_counts") != expected_dispositions:
        errors.append("dataset_profile.disposition_counts does not match records")
    if profile.get("human_review_count") != expected_dispositions[
        "review_required"
    ]:
        errors.append("dataset_profile.human_review_count does not match records")
    statistics = profile.get("descriptor_statistics")
    if not isinstance(statistics, dict) or set(statistics) != descriptor_names:
        errors.append(
            "dataset_profile.descriptor_statistics does not match descriptor_set"
        )
    else:
        for name, stats in statistics.items():
            path = f"dataset_profile.descriptor_statistics.{name}"
            if not isinstance(stats, dict):
                errors.append(f"{path} must be an object")
                continue
            if stats.get("total_records") != len(records):
                errors.append(f"{path}.total_records does not match records")
            missing = stats.get("missing_count")
            non_missing = stats.get("non_missing_count")
            if not isinstance(missing, int) or not isinstance(non_missing, int):
                errors.append(f"{path} missing/non-missing counts must be integers")
            elif missing + non_missing != len(records):
                errors.append(f"{path} counts do not conserve total records")
            expected_rate = missing / len(records) if records else None
            rate = stats.get("missing_rate")
            if expected_rate is not None and (
                not isinstance(rate, (int, float))
                or not math.isclose(
                    float(rate), expected_rate, rel_tol=0.0, abs_tol=1e-15
                )
            ):
                errors.append(f"{path}.missing_rate does not match counts")


def validate(document: Any) -> dict[str, Any]:
    errors: list[str] = []
    warnings: list[str] = []
    if not isinstance(document, dict):
        return {
            "valid": False,
            "errors": ["document must be an object"],
            "warnings": [],
        }
    require_keys(document, REQUIRED_TOP_LEVEL, "document", errors)
    if document.get("schema_version") != SCHEMA_VERSION:
        errors.append(f"schema_version must be {SCHEMA_VERSION}")
    if document.get("workflow") != WORKFLOW:
        errors.append(f"workflow must be {WORKFLOW}")
    if not document.get("generated_at_utc"):
        errors.append("generated_at_utc is required")

    tool_versions = document.get("tool_versions")
    if not isinstance(tool_versions, dict):
        errors.append("tool_versions must be an object")
    else:
        if tool_versions.get("rdkit") not in {"2025.9.2", "2025.09.2"}:
            errors.append("tool_versions.rdkit must be fixed to 2025.9.2")
        for field in ("python", "feature_calculator"):
            if not tool_versions.get(field):
                errors.append(f"tool_versions.{field} is required")

    options = document.get("options")
    if not isinstance(options, dict):
        errors.append("options must be an object")
        calculation_view = "invalid"
    else:
        calculation_view = options.get("calculation_view")
        if calculation_view not in CALCULATION_VIEWS:
            errors.append("options.calculation_view is invalid")
        if options.get("auto_repair_structures") is not False:
            errors.append("options.auto_repair_structures must be false")
        if options.get("requires_3d_conformer") is not False:
            errors.append("options.requires_3d_conformer must be false")

    descriptor_metadata = document.get("descriptor_set")
    descriptor_names: set[str] = set()
    if not isinstance(descriptor_metadata, dict):
        errors.append("descriptor_set must be an object")
    else:
        if descriptor_metadata.get("id") != DESCRIPTOR_SET_ID:
            errors.append("descriptor_set.id is invalid")
        if descriptor_metadata.get("requires_3d_conformer") is not False:
            errors.append("descriptor_set must not require a 3D conformer")
        features = descriptor_metadata.get("features")
        if not isinstance(features, list):
            errors.append("descriptor_set.features must be a list")
        else:
            for index, feature in enumerate(features):
                if not isinstance(feature, dict) or not feature.get("name"):
                    errors.append(
                        f"descriptor_set.features[{index}] is invalid"
                    )
                    continue
                descriptor_names.add(feature["name"])
                if feature.get("feature_class") not in {
                    "structure_deterministic_calculation",
                    "structure_based_empirical_descriptor",
                }:
                    errors.append(
                        f"descriptor_set.features[{index}].feature_class is invalid"
                    )
            if descriptor_names != DESCRIPTOR_NAMES:
                errors.append("descriptor_set features do not match the core set")

    profiles = document.get("fingerprint_profiles")
    if not isinstance(profiles, dict) or set(profiles) != FINGERPRINT_NAMES:
        errors.append("fingerprint_profiles must contain the three core profiles")
        profiles = {}
    else:
        for name, profile in profiles.items():
            validate_profile(
                name, profile, f"fingerprint_profiles.{name}", errors
            )

    records = document.get("records")
    if not isinstance(records, list) or not records:
        errors.append("records must be a non-empty list")
        records = []
    for index, record in enumerate(records):
        validate_record(
            record,
            index,
            descriptor_names,
            profiles,
            calculation_view,
            errors,
            warnings,
        )

    summary = document.get("input_summary")
    if not isinstance(summary, dict):
        errors.append("input_summary must be an object")
    else:
        if summary.get("total_records") != len(records):
            errors.append("input_summary.total_records does not match records")
        expected_statuses = {
            status: sum(
                record.get("calculation_status") == status
                for record in records
                if isinstance(record, dict)
            )
            for status in CALCULATION_STATUSES
        }
        if summary.get("calculation_status_counts") != expected_statuses:
            errors.append(
                "input_summary.calculation_status_counts does not match records"
            )
        expected_dispositions = {
            status: sum(
                record.get("disposition") == status
                for record in records
                if isinstance(record, dict)
            )
            for status in DISPOSITIONS
        }
        if summary.get("output_disposition_counts") != expected_dispositions:
            errors.append(
                "input_summary.output_disposition_counts does not match records"
            )
        if sum(expected_dispositions.values()) != len(records):
            errors.append("record dispositions do not conserve input count")

    validate_dataset_profile(
        document.get("dataset_profile"), records, descriptor_names, errors
    )

    for field in (
        "errors",
        "warnings",
        "notices",
        "human_review_required",
    ):
        if not isinstance(document.get(field), list):
            errors.append(f"{field} must be a list")
    for field in ("errors", "warnings", "human_review_required"):
        for index, item in enumerate(document.get(field, [])):
            validate_finding(item, f"{field}[{index}]", errors)
            if isinstance(item, dict) and not item.get("record_id"):
                errors.append(f"{field}[{index}].record_id is required")

    non_finite_paths = find_non_finite(document)
    if non_finite_paths:
        errors.append(
            "output contains non-finite numeric values at: "
            + ", ".join(non_finite_paths[:20])
        )
    serialized = json.dumps(document, ensure_ascii=False)
    if SECRET_RE.search(serialized):
        errors.append("possible secret detected in output")
    lowered = serialized.lower()
    for claim in FORBIDDEN_CLAIMS:
        if claim.lower() in lowered:
            errors.append(f"forbidden scientific claim detected: {claim}")

    fingerprint = document.get("result_fingerprint")
    if not isinstance(fingerprint, str) or not re.fullmatch(
        r"[0-9a-f]{64}", fingerprint or ""
    ):
        errors.append("result_fingerprint must be a SHA-256 hex string")
    elif fingerprint != output_fingerprint(document):
        errors.append("result_fingerprint mismatch")

    return {
        "valid": not errors,
        "errors": errors,
        "warnings": sorted(set(warnings)),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("path", type=Path)
    args = parser.parse_args()
    try:
        document = json.loads(
            args.path.read_text(encoding="utf-8"),
            parse_constant=lambda value: float(value),
        )
        report = validate(document)
    except (OSError, json.JSONDecodeError, ValueError) as error:
        report = {"valid": False, "errors": [str(error)], "warnings": []}
    sys.stdout.write(json.dumps(report, ensure_ascii=False, indent=2) + "\n")
    return 0 if report["valid"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
