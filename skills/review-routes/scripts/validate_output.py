#!/usr/bin/env python3
"""Validate review-routes output and scientific boundaries."""

from __future__ import annotations

import argparse
import json
import re
from collections import defaultdict
from pathlib import Path
from typing import Any, Sequence

from review_routes import (
    DISPOSITIONS,
    PRECEDENT_LEVELS,
    REVIEW_STATUSES,
    RULE_MESSAGES,
    SCHEMA_VERSION,
    SECRET_RE,
    WORKFLOW,
    stable_document_fingerprint,
)

FORBIDDEN_KEYS = {
    "route_is_feasible",
    "safe_to_execute",
    "optimal_route",
    "recommended_route",
    "success_probability",
    "ready_for_experiment",
    "overall_yield",
    "decision_score",
    "total_score",
}
FORBIDDEN_CLAIMS = {
    "该路线可行",
    "路线安全",
    "最佳路线",
    "可直接实验",
    "route is feasible",
    "safe to execute",
    "optimal route",
    "ready to run",
}
PARTIAL_CODES = {
    "W-CURATION-REVIEW-001",
    "W-PRECEDENT-TIMEOUT-001",
    "W-PRECEDENT-ERROR-001",
    "W-PRECEDENT-NOT-RUN-001",
}


def walk(value: Any, path: str = "$") -> list[tuple[str, str, Any]]:
    output = []
    if isinstance(value, dict):
        for key, item in value.items():
            current = f"{path}.{key}"
            output.append((current, str(key), item))
            output.extend(walk(item, current))
    elif isinstance(value, list):
        for index, item in enumerate(value):
            output.extend(walk(item, f"{path}[{index}]"))
    return output


def validate_finding(value: Any, path: str) -> list[str]:
    errors = []
    if not isinstance(value, dict):
        return [f"{path} 必须是 object"]
    if value.get("code") not in RULE_MESSAGES:
        errors.append(f"{path}.code 未登记")
    if value.get("severity") not in {"error", "warning"}:
        errors.append(f"{path}.severity 不受控")
    if value.get("message") != RULE_MESSAGES.get(value.get("code")):
        errors.append(f"{path}.message 与规则目录不一致")
    if not isinstance(value.get("field_path"), str) or not value["field_path"]:
        errors.append(f"{path}.field_path 为空")
    if not isinstance(value.get("evidence"), list):
        errors.append(f"{path}.evidence 必须是 array")
    return errors


def validate_step(
    value: Any, path: str, *, allow_missing_hash: bool = False
) -> list[str]:
    errors = []
    if not isinstance(value, dict):
        return [f"{path} 必须是 object"]
    required = {
        "step_id",
        "step_reaction_hash",
        "path",
        "reported_reaction",
        "canonical_reaction",
        "product",
        "precursors",
        "agents",
        "backend_metadata",
        "curation",
        "precedent",
        "findings",
        "review_required",
    }
    for key in sorted(required):
        if key not in value:
            errors.append(f"{path}.{key} 缺失")
    if not isinstance(value.get("step_id"), str) or not value.get("step_id"):
        errors.append(f"{path}.step_id 非法")
    if not (
        allow_missing_hash and value.get("step_reaction_hash") is None
    ) and not re.fullmatch(
        r"[0-9a-f]{64}", str(value.get("step_reaction_hash", ""))
    ):
        errors.append(f"{path}.step_reaction_hash 非 SHA-256")
    for key in ("path", "precursors", "agents", "findings", "review_required"):
        if not isinstance(value.get(key), list):
            errors.append(f"{path}.{key} 必须是 array")
    curation = value.get("curation")
    if not isinstance(curation, dict):
        errors.append(f"{path}.curation 必须是 object")
    precedent = value.get("precedent")
    if not isinstance(precedent, dict):
        errors.append(f"{path}.precedent 必须是 object")
    elif precedent.get("match_level") not in PRECEDENT_LEVELS:
        errors.append(f"{path}.precedent.match_level 不受控")
    step_findings = value.get("findings") or []
    for index, item in enumerate(step_findings):
        errors.extend(validate_finding(item, f"{path}.findings[{index}]"))
    expected_codes = sorted(
        {
            item.get("code")
            for item in step_findings
            if isinstance(item, dict) and item.get("code")
        }
    )
    if value.get("review_required") != expected_codes:
        errors.append(f"{path}.review_required 与 findings 不一致")
    return errors


def validate_route(value: Any, path: str) -> list[str]:
    errors = []
    if not isinstance(value, dict):
        return [f"{path} 必须是 object"]
    required = {
        "route_id",
        "source_route_hash",
        "backend_metadata",
        "route_signature",
        "target_structure",
        "topology_status",
        "node_count",
        "step_count",
        "longest_linear_sequence",
        "branch_count",
        "terminal_precursors",
        "inventory_snapshot",
        "inventory_coverage",
        "precedent_coverage_by_level",
        "exact_or_transformation_coverage",
        "weakest_steps",
        "constraint_results",
        "step_reviews",
        "findings",
        "review_status",
        "disposition",
        "human_review_required",
        "duplicate_memberships",
    }
    for key in sorted(required):
        if key not in value:
            errors.append(f"{path}.{key} 缺失")
    if not isinstance(value.get("route_id"), str) or not value.get("route_id"):
        errors.append(f"{path}.route_id 非法")
    if not re.fullmatch(r"[0-9a-f]{64}", str(value.get("source_route_hash", ""))):
        errors.append(f"{path}.source_route_hash 非 SHA-256")
    if not re.fullmatch(r"route:[0-9a-f]{24}", str(value.get("route_signature", ""))):
        errors.append(f"{path}.route_signature 非受控格式")
    if value.get("topology_status") not in {"valid", "invalid"}:
        errors.append(f"{path}.topology_status 不受控")
    if value.get("review_status") not in REVIEW_STATUSES:
        errors.append(f"{path}.review_status 不受控")
    if value.get("disposition") not in DISPOSITIONS:
        errors.append(f"{path}.disposition 不受控")
    for key in (
        "terminal_precursors",
        "weakest_steps",
        "constraint_results",
        "step_reviews",
        "findings",
        "human_review_required",
        "duplicate_memberships",
    ):
        if not isinstance(value.get(key), list):
            errors.append(f"{path}.{key} 必须是 array")
    steps = value.get("step_reviews") or []
    if value.get("step_count") != len(steps):
        errors.append(f"{path}.step_count 与 step_reviews 不一致")
    route_findings = value.get("findings") or []
    allow_missing_hash = any(
        isinstance(item, dict)
        and item.get("code") == "E-STEP-REACTION-001"
        and item.get("severity") == "error"
        for item in route_findings
    )
    step_ids = []
    for index, item in enumerate(steps):
        errors.extend(
            validate_step(
                item,
                f"{path}.step_reviews[{index}]",
                allow_missing_hash=allow_missing_hash,
            )
        )
        if isinstance(item, dict):
            step_ids.append(item.get("step_id"))
    if len(step_ids) != len(set(step_ids)):
        errors.append(f"{path}.step_id 批内重复")
    coverage = value.get("precedent_coverage_by_level")
    if not isinstance(coverage, dict):
        errors.append(f"{path}.precedent_coverage_by_level 必须是 object")
    else:
        if set(coverage) != PRECEDENT_LEVELS:
            errors.append(f"{path}.precedent_coverage_by_level 枚举不完整")
        if any(not isinstance(item, int) or item < 0 for item in coverage.values()):
            errors.append(f"{path}.precedent_coverage_by_level 必须是非负整数")
        if sum(coverage.values()) != len(steps):
            errors.append(f"{path}.precedent coverage 计数不守恒")
    exact_coverage = value.get("exact_or_transformation_coverage")
    if (
        not isinstance(exact_coverage, (int, float))
        or isinstance(exact_coverage, bool)
        or not 0 <= exact_coverage <= 1
    ):
        errors.append(f"{path}.exact_or_transformation_coverage 非 0–1")
    inventory_coverage = value.get("inventory_coverage")
    if (
        not isinstance(inventory_coverage, (int, float))
        or isinstance(inventory_coverage, bool)
        or not 0 <= inventory_coverage <= 1
    ):
        errors.append(f"{path}.inventory_coverage 非 0–1")
    findings = value.get("findings") or []
    for index, item in enumerate(findings):
        errors.extend(validate_finding(item, f"{path}.findings[{index}]"))
    severities = {
        item.get("severity") for item in findings if isinstance(item, dict)
    }
    codes = sorted(
        {
            item.get("code")
            for item in findings
            if isinstance(item, dict) and item.get("code")
        }
    )
    expected_status = (
        "error"
        if "error" in severities
        else "partial"
        if any(code in PARTIAL_CODES for code in codes)
        else "completed"
    )
    expected_disposition = (
        "blocked"
        if "error" in severities
        else "review_required"
        if findings
        else "ready_for_expert_review"
    )
    if value.get("review_status") != expected_status:
        errors.append(f"{path}.review_status 应为 {expected_status}")
    if value.get("disposition") != expected_disposition:
        errors.append(f"{path}.disposition 应为 {expected_disposition}")
    expected_human = [code for code in codes if code.startswith("W-")]
    if value.get("human_review_required") != expected_human:
        errors.append(f"{path}.human_review_required 与 findings 不一致")
    if (
        any(item.get("code") == "E-CURATION-REJECTED-001" for item in findings)
        and value.get("disposition") != "blocked"
    ):
        errors.append(f"{path} rejected step 不得进入 ready/review")
    return errors


def validate_output(document: Any) -> list[str]:
    errors = []
    if not isinstance(document, dict):
        return ["输出顶层必须是 object"]
    required = {
        "schema_version",
        "workflow",
        "ruleset_version",
        "generated_at_utc",
        "tool_versions",
        "source_record",
        "routes_fingerprint",
        "options",
        "constraints",
        "target_assessment",
        "input_summary",
        "route_summaries",
        "duplicate_route_groups",
        "comparison_dimensions",
        "review_queue",
        "errors",
        "warnings",
        "notices",
        "result_fingerprint",
    }
    for key in sorted(required):
        if key not in document:
            errors.append(f"{key} 缺失")
    if document.get("schema_version") != SCHEMA_VERSION:
        errors.append("schema_version 不匹配")
    if document.get("workflow") != WORKFLOW:
        errors.append("workflow 不匹配")
    versions = document.get("tool_versions")
    if not isinstance(versions, dict) or versions.get("rdkit") not in {
        "2025.9.2",
        "2025.09.2",
    }:
        errors.append("rdkit 必须固定 2025.9.2")
    options = document.get("options")
    expected_options = {
        "comparison_mode": "dimensions_only",
        "preserve_backend_order": True,
        "automatic_route_ranking": False,
        "network_access": False,
        "pickle_allowed": False,
    }
    if not isinstance(options, dict):
        errors.append("options 必须是 object")
    else:
        for key, value in expected_options.items():
            if options.get(key) != value:
                errors.append(f"options.{key} 必须是 {value!r}")
    for key in (
        "route_summaries",
        "duplicate_route_groups",
        "comparison_dimensions",
        "review_queue",
        "errors",
        "warnings",
        "notices",
    ):
        if not isinstance(document.get(key), list):
            errors.append(f"{key} 必须是 array")
    routes = document.get("route_summaries") or []
    route_ids = []
    for index, item in enumerate(routes):
        errors.extend(validate_route(item, f"route_summaries[{index}]"))
        if isinstance(item, dict):
            route_ids.append(item.get("route_id"))
    if len(route_ids) != len(set(route_ids)):
        errors.append("route_id 输出重复")
    summary = document.get("input_summary")
    if not isinstance(summary, dict):
        errors.append("input_summary 必须是 object")
    else:
        if summary.get("output_routes") != len(routes):
            errors.append("input_summary.output_routes 不一致")
        if summary.get("record_count_conserved") != (
            summary.get("input_routes") == summary.get("output_routes")
        ):
            errors.append("input_summary.record_count_conserved 不一致")
        counts = summary.get("disposition_counts")
        if not isinstance(counts, dict) or set(counts) != DISPOSITIONS:
            errors.append("input_summary.disposition_counts 不完整")
        elif sum(counts.values()) != len(routes):
            errors.append("input_summary.disposition_counts 不守恒")
    comparisons = document.get("comparison_dimensions") or []
    if len(comparisons) != len(routes):
        errors.append("comparison_dimensions 与 routes 数量不一致")
    comparison_ids = {
        item.get("route_id") for item in comparisons if isinstance(item, dict)
    }
    if comparison_ids != set(route_ids):
        errors.append("comparison_dimensions route_id 不一致")
    groups = document.get("duplicate_route_groups") or []
    memberships = defaultdict(set)
    for index, group in enumerate(groups):
        if not isinstance(group, dict):
            errors.append(f"duplicate_route_groups[{index}] 必须是 object")
            continue
        group_id = group.get("group_id")
        members = group.get("route_ids")
        if not isinstance(members, list) or len(members) < 2:
            errors.append(f"duplicate_route_groups[{index}].route_ids 非重复组")
            continue
        for route_id in members:
            memberships[route_id].add(group_id)
    for route in routes:
        if isinstance(route, dict) and set(route.get("duplicate_memberships") or []) != memberships.get(
            route.get("route_id"), set()
        ):
            errors.append(
                f"route {route.get('route_id')} duplicate_memberships 不一致"
            )
    for path, key, _ in walk(document):
        if key in FORBIDDEN_KEYS:
            errors.append(f"{path} 是禁止字段")
    serialized = json.dumps(document, ensure_ascii=False)
    if SECRET_RE.search(serialized):
        errors.append("输出含疑似凭证")
    for claim in FORBIDDEN_CLAIMS:
        if re.search(re.escape(claim), serialized, flags=re.IGNORECASE):
            errors.append(f"输出含禁止科学结论：{claim}")
    if document.get("result_fingerprint") != stable_document_fingerprint(document):
        errors.append("result_fingerprint 不匹配")
    return errors


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    args = parser.parse_args(argv)
    try:
        document = json.loads(args.output.read_text(encoding="utf-8"))
    except Exception as error:
        print(json.dumps({"valid": False, "errors": [str(error)]}, ensure_ascii=False))
        return 2
    errors = validate_output(document)
    print(
        json.dumps(
            {"valid": not errors, "errors": errors},
            ensure_ascii=False,
            indent=2,
        )
    )
    return 0 if not errors else 1


if __name__ == "__main__":
    raise SystemExit(main())
