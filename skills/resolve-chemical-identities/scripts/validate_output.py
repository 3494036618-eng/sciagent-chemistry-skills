#!/usr/bin/env python3
"""校验 resolve-chemical-identities 的输出契约和科学失败关闭规则。"""

from __future__ import annotations

import argparse
import importlib.util
import json
import re
import sys
from pathlib import Path
from typing import Any


SCHEMA_VERSION = "1.0.0"
WORKFLOW = "chemical-identity-resolution"
INPUT_STATUSES = {"valid", "invalid_input"}
RETRIEVAL_STATUSES = {
    "completed",
    "partial",
    "not_found",
    "source_error",
    "not_run",
}
ALIGNMENT_STATUSES = {
    "exact",
    "related_forms",
    "ambiguous",
    "conflict",
    "not_assessed",
}
SAMPLE_STATUSES = {"not_assessed", "user_confirmed", "expert_confirmed"}
DISPOSITIONS = {
    "ready_for_standardization",
    "review_required",
    "rejected",
}
SOURCE_QUERY_STATUSES = {"success", "not_found", "source_error"}
SECRET_RE = re.compile(
    r"(?i)(authorization\s*:|bearer\s+[A-Za-z0-9._~+/=-]{12,}|"
    r"api[_ -]?key\s*[:=]|cookie\s*:|ark-[A-Za-z0-9_-]{16,})"
)
FORBIDDEN_CLAIMS = (
    "sample identity confirmed",
    "physical sample confirmed",
    "experimentally confirmed",
    "safe to synthesize",
    "proven efficacy",
    "proven active",
)


def load_resolver() -> Any:
    path = Path(__file__).with_name("resolve_identities.py")
    spec = importlib.util.spec_from_file_location("_identity_resolver", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"无法加载 resolver：{path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def require_keys(
    value: dict[str, Any],
    keys: set[str],
    path: str,
    errors: list[str],
) -> None:
    missing = sorted(keys - set(value))
    if missing:
        errors.append(f"{path} missing keys: {', '.join(missing)}")


def valid_enum(
    value: Any, allowed: set[str], path: str, errors: list[str]
) -> None:
    if value not in allowed:
        errors.append(f"{path} has invalid value: {value!r}")


def validate_candidate(
    candidate: Any,
    path: str,
    errors: list[str],
    warnings: list[str],
) -> None:
    if not isinstance(candidate, dict):
        errors.append(f"{path} must be an object")
        return
    require_keys(
        candidate,
        {
            "candidate_id",
            "canonical_smiles",
            "inchi",
            "inchikey",
            "connectivity_block",
            "names",
            "evidence",
            "source_families",
            "quality_findings",
            "comparison_view",
        },
        path,
        errors,
    )
    inchikey = candidate.get("inchikey")
    if not isinstance(inchikey, str) or len(inchikey) != 27:
        errors.append(f"{path}.inchikey must be a full InChIKey")
    if not isinstance(candidate.get("evidence"), list) or not candidate.get("evidence"):
        errors.append(f"{path}.evidence must be a non-empty list")
    if not isinstance(candidate.get("source_families"), list):
        errors.append(f"{path}.source_families must be a list")
    view = candidate.get("comparison_view")
    if view is not None and not isinstance(view, dict):
        errors.append(f"{path}.comparison_view must be null or an object")
    if isinstance(view, dict):
        if view.get("parent_inchikey") and not view.get("parent_structure"):
            errors.append(
                f"{path}.comparison_view parent_inchikey requires parent_structure"
            )
        if view.get("status") == "completed" and not view.get(
            "standardized_structure"
        ):
            errors.append(
                f"{path}.comparison_view completed without standardized_structure"
            )
    if not candidate.get("canonical_smiles"):
        warnings.append(f"{path} has no locally reconstructed canonical SMILES")


def validate_resolution(
    resolution: Any,
    index: int,
    errors: list[str],
    warnings: list[str],
) -> None:
    path = f"resolutions[{index}]"
    if not isinstance(resolution, dict):
        errors.append(f"{path} must be an object")
        return
    require_keys(
        resolution,
        {
            "request",
            "input_status",
            "retrieval_status",
            "record_alignment_status",
            "record_alignment_scope",
            "sample_identity_status",
            "disposition",
            "candidates",
            "unresolved_source_records",
            "source_record_conflicts",
            "source_queries",
            "relationship_evidence",
            "confirmation_questions",
            "findings",
            "standardization_comparison",
            "standardization_handoff",
        },
        path,
        errors,
    )
    valid_enum(resolution.get("input_status"), INPUT_STATUSES, f"{path}.input_status", errors)
    valid_enum(
        resolution.get("retrieval_status"),
        RETRIEVAL_STATUSES,
        f"{path}.retrieval_status",
        errors,
    )
    valid_enum(
        resolution.get("record_alignment_status"),
        ALIGNMENT_STATUSES,
        f"{path}.record_alignment_status",
        errors,
    )
    valid_enum(
        resolution.get("sample_identity_status"),
        SAMPLE_STATUSES,
        f"{path}.sample_identity_status",
        errors,
    )
    valid_enum(
        resolution.get("disposition"),
        DISPOSITIONS,
        f"{path}.disposition",
        errors,
    )

    request = resolution.get("request")
    if not isinstance(request, dict) or not isinstance(request.get("query"), str):
        errors.append(f"{path}.request.query must preserve the original string")
    if resolution.get("record_alignment_scope") != "database_records_only":
        errors.append(f"{path}.record_alignment_scope must be database_records_only")
    if resolution.get("sample_identity_status") != "not_assessed":
        errors.append(
            f"{path}.sample_identity_status cannot be automatically upgraded"
        )

    candidates = resolution.get("candidates")
    if not isinstance(candidates, list):
        errors.append(f"{path}.candidates must be a list")
        candidates = []
    for candidate_index, candidate in enumerate(candidates):
        validate_candidate(
            candidate,
            f"{path}.candidates[{candidate_index}]",
            errors,
            warnings,
        )

    source_queries = resolution.get("source_queries")
    if not isinstance(source_queries, list):
        errors.append(f"{path}.source_queries must be a list")
        source_queries = []
    for query_index, query in enumerate(source_queries):
        query_path = f"{path}.source_queries[{query_index}]"
        if not isinstance(query, dict):
            errors.append(f"{query_path} must be an object")
            continue
        valid_enum(
            query.get("status"),
            SOURCE_QUERY_STATUSES,
            f"{query_path}.status",
            errors,
        )
        if query.get("status") == "source_error" and not query.get("error_kind"):
            errors.append(f"{query_path} source_error requires error_kind")
        if query.get("status") == "not_found" and query.get("http_status") not in {
            None,
            200,
            404,
        }:
            warnings.append(
                f"{query_path} not_found did not originate from HTTP 404 or an empty 200 response"
            )

    statuses = [item.get("status") for item in source_queries if isinstance(item, dict)]
    retrieval = resolution.get("retrieval_status")
    if retrieval == "not_found" and any(status == "source_error" for status in statuses):
        errors.append(f"{path} misclassifies source_error as not_found")
    if retrieval == "source_error" and any(status == "success" for status in statuses):
        errors.append(f"{path} must use partial when success and source_error coexist")
    if retrieval == "partial" and not (
        "success" in statuses and "source_error" in statuses
    ):
        errors.append(f"{path}.partial requires both success and source_error")

    alignment = resolution.get("record_alignment_status")
    disposition = resolution.get("disposition")
    if alignment == "exact" and len(candidates) != 1:
        errors.append(f"{path}.exact requires exactly one complete candidate")
    if alignment in {"ambiguous", "conflict", "related_forms"}:
        if disposition != "review_required":
            errors.append(f"{path}.{alignment} must require review")
        if not resolution.get("confirmation_questions"):
            errors.append(f"{path}.{alignment} requires a confirmation question")
    if alignment == "exact" and candidates:
        input_type = (request or {}).get("detected_input_type")
        families = set(candidates[0].get("source_families") or [])
        if input_type == "name":
            independent_external = families - {"local_input", "unichem"}
            if len(independent_external) < 2:
                errors.append(
                    f"{path}.exact name requires at least two independent external sources"
                )
    if resolution.get("source_record_conflicts") and alignment != "conflict":
        errors.append(f"{path} source record conflict must set alignment=conflict")

    if resolution.get("input_status") == "invalid_input":
        if retrieval != "not_run" or candidates:
            errors.append(f"{path} invalid input must not query sources or emit candidates")
        if disposition != "rejected":
            errors.append(f"{path} invalid input must be rejected")

    handoff = resolution.get("standardization_handoff")
    if not isinstance(handoff, dict):
        errors.append(f"{path}.standardization_handoff must be an object")
    elif handoff.get("status") == "ready":
        if disposition != "ready_for_standardization" or len(candidates) != 1:
            errors.append(f"{path} ready handoff requires one ready candidate")
        if not handoff.get("records"):
            errors.append(f"{path} ready handoff requires records")
    elif handoff.get("records"):
        errors.append(f"{path} blocked handoff must not contain records")


def validate(document: Any) -> dict[str, Any]:
    errors: list[str] = []
    warnings: list[str] = []
    if not isinstance(document, dict):
        return {
            "valid": False,
            "errors": ["document must be an object"],
            "warnings": [],
        }
    require_keys(
        document,
        {
            "schema_version",
            "workflow",
            "generated_at_utc",
            "tool_versions",
            "source_metadata",
            "options",
            "input_summary",
            "resolutions",
            "cross_query_relationships",
            "notices",
            "result_fingerprint",
        },
        "document",
        errors,
    )
    if document.get("schema_version") != SCHEMA_VERSION:
        errors.append("unsupported schema_version")
    if document.get("workflow") != WORKFLOW:
        errors.append("invalid workflow")
    source_metadata = document.get("source_metadata")
    if not isinstance(source_metadata, dict):
        errors.append("source_metadata must be an object")
    else:
        for source in ("OPSIN", "PubChem", "ChEMBL", "UniChem"):
            if source not in source_metadata:
                errors.append(f"source_metadata missing {source}")
    options = document.get("options")
    if not isinstance(options, dict):
        errors.append("options must be an object")
    else:
        if options.get("automatic_tie_breaking") is not False:
            errors.append("automatic_tie_breaking must be false")
        if options.get("no_model_generated_structures") is not True:
            errors.append("no_model_generated_structures must be true")

    resolutions = document.get("resolutions")
    if not isinstance(resolutions, list) or not resolutions:
        errors.append("resolutions must be a non-empty list")
        resolutions = []
    for index, resolution in enumerate(resolutions):
        validate_resolution(resolution, index, errors, warnings)

    summary = document.get("input_summary")
    if not isinstance(summary, dict):
        errors.append("input_summary must be an object")
    else:
        counts = {
            status: sum(
                resolution.get("disposition") == status
                for resolution in resolutions
                if isinstance(resolution, dict)
            )
            for status in DISPOSITIONS
        }
        if summary.get("total_requests") != len(resolutions):
            errors.append("input_summary.total_requests does not match resolutions")
        for status, count in counts.items():
            if summary.get(status) != count:
                errors.append(f"input_summary.{status} does not match resolutions")
        if sum(counts.values()) != len(resolutions):
            errors.append("disposition counts do not conserve total requests")

    serialized = json.dumps(document, ensure_ascii=False)
    if SECRET_RE.search(serialized):
        errors.append("possible secret detected in output")
    lower = serialized.lower()
    for claim in FORBIDDEN_CLAIMS:
        if claim in lower:
            errors.append(f"forbidden scientific claim detected: {claim}")

    try:
        resolver = load_resolver()
        expected = resolver.output_fingerprint(document)
        if document.get("result_fingerprint") != expected:
            errors.append("result_fingerprint mismatch")
    except Exception as error:
        errors.append(f"could not validate result_fingerprint: {error}")

    return {
        "valid": not errors,
        "errors": errors,
        "warnings": warnings,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("path", type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        document = json.loads(args.path.read_text(encoding="utf-8"))
        report = validate(document)
    except (OSError, json.JSONDecodeError) as error:
        report = {"valid": False, "errors": [str(error)], "warnings": []}
    sys.stdout.write(json.dumps(report, ensure_ascii=False, indent=2) + "\n")
    return 0 if report["valid"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
