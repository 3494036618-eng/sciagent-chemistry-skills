#!/usr/bin/env python3
"""校验 standardize-chemical-structures 的 JSON 输出契约。"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path
from typing import Any


SCHEMA_VERSION = "1.0.0"
WORKFLOW = "chemical-structure-standardization-qc"
PROFILES = {"rdkit-basic", "chembl-pipeline"}
DISPOSITIONS = {"ready_for_downstream", "review_required", "rejected"}
PARSE_STATUSES = {"success", "error"}
STANDARDIZATION_STATUSES = {"completed", "not_run", "error"}
DUPLICATE_BASES = {"original", "standardized", "parent"}
REQUIRED_TOP_LEVEL = {
    "schema_version",
    "workflow",
    "generated_at_utc",
    "tool_versions",
    "options",
    "input_summary",
    "records",
    "duplicate_groups",
    "errors",
    "warnings",
    "notices",
    "human_review_required",
    "result_fingerprint",
}
REQUIRED_RECORD_FIELDS = {
    "id",
    "original_structure",
    "input_format",
    "parse_status",
    "standardization_status",
    "standardized_structure",
    "parent_structure",
    "inchikey",
    "parent_inchikey",
    "transformations",
    "qc_findings",
    "disposition",
    "human_review_required",
}
SECRET_RE = re.compile(
    r"ark-[A-Za-z0-9_-]{12,}|"
    r"Bearer\s+[A-Za-z0-9._-]{12,}|"
    r"(?:Authorization|Cookie|Token)\s*[:=]\s*[A-Za-z0-9._-]{12,}",
    re.IGNORECASE,
)
FORBIDDEN_ASSERTIONS = {
    "实验样品身份已经确认",
    "结构已确证",
    "药效已确认",
    "毒性已确认",
    "安全性已确认",
    "可合成性已确认",
    "experimentally confirmed",
    "clinically effective",
    "safe to synthesize",
}


def output_fingerprint(document: dict[str, Any]) -> str:
    payload = {
        key: value
        for key, value in document.items()
        if key not in {"generated_at_utc", "result_fingerprint"}
    }
    serialized = json.dumps(
        payload, ensure_ascii=False, sort_keys=True, separators=(",", ":")
    )
    return hashlib.sha256(serialized.encode("utf-8")).hexdigest()


def validate_finding(
    item: Any,
    path: str,
    errors: list[str],
) -> None:
    if not isinstance(item, dict):
        errors.append(f"{path} must be an object")
        return
    for field in ("code", "severity", "message", "source"):
        if not item.get(field):
            errors.append(f"{path}.{field} is required")
    if item.get("severity") not in {"error", "warning", "review"}:
        errors.append(f"{path}.severity is invalid")


def validate_record(
    record: Any,
    index: int,
    errors: list[str],
    warnings: list[str],
) -> None:
    path = f"records[{index}]"
    if not isinstance(record, dict):
        errors.append(f"{path} must be an object")
        return
    missing = REQUIRED_RECORD_FIELDS - set(record)
    if missing:
        errors.append(f"{path} missing fields: {sorted(missing)!r}")
        return
    if not isinstance(record["id"], str) or not record["id"]:
        errors.append(f"{path}.id must be a non-empty string")
    if not isinstance(record["original_structure"], str):
        errors.append(f"{path}.original_structure must be a string")
    if record["input_format"] not in {"smiles", "sdf", "molblock"}:
        errors.append(f"{path}.input_format is invalid")
    if record["parse_status"] not in PARSE_STATUSES:
        errors.append(f"{path}.parse_status is invalid")
    if record["standardization_status"] not in STANDARDIZATION_STATUSES:
        errors.append(f"{path}.standardization_status is invalid")
    if record["disposition"] not in DISPOSITIONS:
        errors.append(f"{path}.disposition is invalid")
    if not isinstance(record["transformations"], list):
        errors.append(f"{path}.transformations must be a list")
    if not isinstance(record["qc_findings"], list):
        errors.append(f"{path}.qc_findings must be a list")
    else:
        for finding_index, item in enumerate(record["qc_findings"]):
            validate_finding(
                item,
                f"{path}.qc_findings[{finding_index}]",
                errors,
            )
    if not isinstance(record["human_review_required"], list):
        errors.append(f"{path}.human_review_required must be a list")

    if record["parse_status"] == "error":
        for field in (
            "standardized_structure",
            "parent_structure",
            "inchikey",
            "parent_inchikey",
        ):
            if record[field] is not None:
                errors.append(
                    f"{path}.{field} must be null when parsing failed"
                )
        if record["disposition"] != "rejected":
            errors.append(f"{path} parse failure must be rejected")

    if record["disposition"] == "ready_for_downstream":
        if record["parse_status"] != "success":
            errors.append(f"{path} ready record must parse successfully")
        if record["standardization_status"] != "completed":
            errors.append(f"{path} ready record must complete standardization")
        if record["human_review_required"]:
            errors.append(f"{path} ready record cannot require human review")

    if record["parent_structure"] is None and record["parent_inchikey"] is not None:
        errors.append(f"{path}.parent_inchikey requires parent_structure")

    review_codes = {
        item.get("code")
        for item in record.get("qc_findings", [])
        if isinstance(item, dict) and item.get("severity") == "review"
    }
    if set(record.get("human_review_required", [])) != review_codes:
        errors.append(f"{path}.human_review_required does not match findings")
    if review_codes and record["disposition"] == "ready_for_downstream":
        errors.append(f"{path} review finding cannot be ready")

    chembl_excluded = any(
        isinstance(item, dict)
        and item.get("step") == "chembl_get_parent"
        and item.get("exclusion_flag") is True
        for item in record.get("transformations", [])
    )
    if chembl_excluded and "R-CHEMBL-EXCLUDED" not in review_codes:
        errors.append(
            f"{path} ChEMBL exclusion flag must require human review"
        )

    fragment_analysis = record.get("fragment_analysis")
    if isinstance(fragment_analysis, dict):
        if (
            fragment_analysis.get("classification") == "mixture_or_complex"
            and record["parent_structure"] is not None
        ):
            errors.append(f"{path} mixture must not be collapsed to one parent")

    if record["disposition"] != "ready_for_downstream":
        warnings.append(f"{path} disposition is {record['disposition']}")


def validate_duplicate_groups(
    groups: Any,
    records: list[dict[str, Any]],
    errors: list[str],
) -> None:
    if not isinstance(groups, list):
        errors.append("duplicate_groups must be a list")
        return
    valid_indices = set(range(len(records)))
    for index, group in enumerate(groups):
        path = f"duplicate_groups[{index}]"
        if not isinstance(group, dict):
            errors.append(f"{path} must be an object")
            continue
        if group.get("basis") not in DUPLICATE_BASES:
            errors.append(f"{path}.basis is invalid")
        indices = group.get("record_indices")
        ids = group.get("record_ids")
        if not isinstance(indices, list) or len(indices) < 2:
            errors.append(f"{path}.record_indices must contain at least two")
            continue
        if any(item not in valid_indices for item in indices):
            errors.append(f"{path} references an unknown record index")
        if not isinstance(ids, list) or len(ids) != len(indices):
            errors.append(f"{path}.record_ids must align with indices")
        if group.get("basis") == "parent" and group.get("relationship") != (
            "same_derived_parent_not_same_physical_sample"
        ):
            errors.append(f"{path} must preserve the parent/sample distinction")


def validate(document: dict[str, Any]) -> dict[str, Any]:
    errors: list[str] = []
    warnings: list[str] = []
    if not isinstance(document, dict):
        return {
            "valid": False,
            "errors": ["document must be an object"],
            "warnings": [],
        }

    missing = REQUIRED_TOP_LEVEL - set(document)
    if missing:
        errors.append(f"missing top-level fields: {sorted(missing)!r}")
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
        for field in ("python", "rdkit", "chembl_structure_pipeline"):
            if not tool_versions.get(field):
                errors.append(f"tool_versions.{field} is required")

    options = document.get("options")
    if not isinstance(options, dict):
        errors.append("options must be an object")
    elif options.get("profile") not in PROFILES:
        errors.append("options.profile is invalid")

    records = document.get("records")
    if not isinstance(records, list):
        errors.append("records must be a list")
        records = []
    elif not records:
        errors.append("records must contain at least one input record")
    for index, record in enumerate(records):
        validate_record(record, index, errors, warnings)

    summary = document.get("input_summary")
    if not isinstance(summary, dict):
        errors.append("input_summary must be an object")
    else:
        expected = {
            "ready_for_downstream": 0,
            "review_required": 0,
            "rejected": 0,
        }
        for record in records:
            disposition = record.get("disposition")
            if disposition in expected:
                expected[disposition] += 1
        if summary.get("total_records") != len(records):
            errors.append("input_summary.total_records does not match records")
        for field, count in expected.items():
            if summary.get(field) != count:
                errors.append(f"input_summary.{field} does not match records")
        if len(records) != sum(expected.values()):
            errors.append("record dispositions do not conserve input count")

    validate_duplicate_groups(
        document.get("duplicate_groups"),
        records,
        errors,
    )

    for field in ("errors", "warnings", "notices", "human_review_required"):
        if not isinstance(document.get(field), list):
            errors.append(f"{field} must be a list")
    for field in ("errors", "warnings", "human_review_required"):
        for index, item in enumerate(document.get(field, [])):
            validate_finding(item, f"{field}[{index}]", errors)
            if isinstance(item, dict) and not item.get("record_id"):
                errors.append(f"{field}[{index}].record_id is required")

    serialized = json.dumps(document, ensure_ascii=False)
    if SECRET_RE.search(serialized):
        errors.append("possible secret detected in output")
    lowered = serialized.lower()
    for phrase in FORBIDDEN_ASSERTIONS:
        if phrase.lower() in lowered:
            errors.append(f"forbidden scientific assertion detected: {phrase}")

    fingerprint = document.get("result_fingerprint")
    if not isinstance(fingerprint, str) or not re.fullmatch(
        r"[0-9a-f]{64}", fingerprint or ""
    ):
        errors.append("result_fingerprint must be a SHA-256 hex string")
    elif fingerprint != output_fingerprint(document):
        errors.append("result_fingerprint does not match document content")

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
        document = json.loads(args.path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        print(
            json.dumps(
                {"valid": False, "errors": [str(error)], "warnings": []},
                ensure_ascii=False,
                indent=2,
            )
        )
        return 2
    report = validate(document)
    print(json.dumps(report, ensure_ascii=False, indent=2))
    return 0 if report["valid"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
