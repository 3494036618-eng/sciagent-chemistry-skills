#!/usr/bin/env python3
"""Validate curate-reactions output contracts and scientific boundaries."""

from __future__ import annotations

import argparse
import json
import re
import sys
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from curate_reactions import (
    CURATION_STATUSES,
    DISPOSITIONS,
    RULE_MESSAGES,
    SCHEMA_VERSION,
    SECRET_RE,
    WORKFLOW,
    stable_document_fingerprint,
)

FORBIDDEN_KEYS = {
    "ready_for_modeling",
    "scientifically_correct",
    "safe_to_execute",
    "experiment_is_reproducible",
    "automatic_deletion",
    "automatic_merge",
}
FORBIDDEN_CLAIMS = {
    "该结果适合建模",
    "该反应正确",
    "安全性已确认",
    "可直接执行",
    "保证可复现",
    "ready for modeling",
    "scientifically correct",
    "proven safe",
}


def walk_keys(value: Any, path: str = "$") -> list[tuple[str, str]]:
    results = []
    if isinstance(value, dict):
        for key, item in value.items():
            current = f"{path}.{key}"
            results.append((current, str(key)))
            results.extend(walk_keys(item, current))
    elif isinstance(value, list):
        for index, item in enumerate(value):
            results.extend(walk_keys(item, f"{path}[{index}]"))
    return results


def validate_finding(
    item: Any, path: str, record_id: Any = None
) -> list[str]:
    errors = []
    if not isinstance(item, dict):
        return [f"{path} 必须是 object"]
    code = item.get("code")
    if code not in RULE_MESSAGES:
        errors.append(f"{path}.code 未登记：{code}")
    if item.get("severity") not in {"error", "warning", "human_review"}:
        errors.append(f"{path}.severity 不受控")
    if not isinstance(item.get("field_path"), str) or not item["field_path"]:
        errors.append(f"{path}.field_path 不得为空")
    if item.get("message") != RULE_MESSAGES.get(code):
        errors.append(f"{path}.message 与规则目录不一致")
    if not isinstance(item.get("evidence"), list):
        errors.append(f"{path}.evidence 必须是 array")
    if record_id is not None and item.get("record_id") not in {None, record_id}:
        errors.append(f"{path}.record_id 不一致")
    return errors


def validate_record(record: Any, path: str) -> list[str]:
    errors = []
    if not isinstance(record, dict):
        return [f"{path} 必须是 object"]
    required = {
        "record_id",
        "source_locator",
        "original_record_hash",
        "ord_record",
        "reaction_smiles",
        "participant_assessments",
        "role_assessment",
        "yield_assessment",
        "balance_assessment",
        "mapping_assessment",
        "duplicate_memberships",
        "curation_status",
        "findings",
        "disposition",
        "human_review_required",
    }
    for key in sorted(required):
        if key not in record:
            errors.append(f"{path}.{key} 缺失")
    if record.get("curation_status") not in CURATION_STATUSES:
        errors.append(f"{path}.curation_status 不受控")
    if record.get("disposition") not in DISPOSITIONS:
        errors.append(f"{path}.disposition 不受控")
    if not re.fullmatch(r"[0-9a-f]{64}", str(record.get("original_record_hash", ""))):
        errors.append(f"{path}.original_record_hash 非 SHA-256")
    for key in (
        "participant_assessments",
        "duplicate_memberships",
        "findings",
        "human_review_required",
    ):
        if not isinstance(record.get(key), list):
            errors.append(f"{path}.{key} 必须是 array")
    findings = record.get("findings") or []
    for index, item in enumerate(findings):
        errors.extend(
            validate_finding(
                item, f"{path}.findings[{index}]", record.get("record_id")
            )
        )
    severities = {
        item.get("severity") for item in findings if isinstance(item, dict)
    }
    expected_disposition = (
        "rejected"
        if "error" in severities
        else "review_required"
        if findings
        else "ready_for_search"
    )
    if record.get("disposition") != expected_disposition:
        errors.append(
            f"{path}.disposition 应为 {expected_disposition}"
        )
    expected_status = (
        "error"
        if "error" in severities
        else "partial"
        if findings
        else "completed"
    )
    if record.get("curation_status") != expected_status:
        errors.append(f"{path}.curation_status 应为 {expected_status}")
    human_codes = sorted(
        {
            item.get("code")
            for item in findings
            if isinstance(item, dict) and item.get("severity") == "human_review"
        }
    )
    if record.get("human_review_required") != human_codes:
        errors.append(f"{path}.human_review_required 与 findings 不一致")
    mapping = record.get("mapping_assessment")
    if isinstance(mapping, dict) and mapping.get("status") == "completed":
        errors.append(f"{path}.mapping_assessment 首版不得声称 completed")
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
        "options",
        "source_record",
        "upstream_artifacts",
        "input_summary",
        "records",
        "duplicate_groups",
        "review_queue",
        "errors",
        "warnings",
        "notices",
        "human_review_required",
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
    if not isinstance(versions, dict):
        errors.append("tool_versions 必须是 object")
    else:
        if versions.get("rdkit") != "2025.9.2":
            errors.append("rdkit 必须固定 2025.9.2")
        if versions.get("ord-schema") != "0.8.3":
            errors.append("ord-schema 必须固定 0.8.3")
    options = document.get("options")
    if not isinstance(options, dict):
        errors.append("options 必须是 object")
    else:
        expected = {
            "participant_view": "reported_form",
            "atom_mapping": "off",
            "balance_check": "diagnostic",
            "preserve_original": True,
            "network_access": False,
            "automatic_writeback": False,
        }
        for key, value in expected.items():
            if options.get(key) != value:
                errors.append(f"options.{key} 必须是 {value!r}")
    arrays = (
        "upstream_artifacts",
        "records",
        "duplicate_groups",
        "review_queue",
        "errors",
        "warnings",
        "notices",
        "human_review_required",
    )
    for key in arrays:
        if not isinstance(document.get(key), list):
            errors.append(f"{key} 必须是 array")
    records = document.get("records") or []
    ids = []
    for index, record in enumerate(records):
        errors.extend(validate_record(record, f"records[{index}]"))
        if isinstance(record, dict):
            ids.append(record.get("record_id"))
    if len(ids) != len(set(ids)):
        errors.append("输出 record_id 不唯一")
    summary = document.get("input_summary")
    if not isinstance(summary, dict):
        errors.append("input_summary 必须是 object")
    else:
        if summary.get("output_records") != len(records):
            errors.append("input_summary.output_records 不守恒")
        dispositions = {
            status: sum(
                isinstance(record, dict)
                and record.get("disposition") == status
                for record in records
            )
            for status in sorted(DISPOSITIONS)
        }
        statuses = {
            status: sum(
                isinstance(record, dict)
                and record.get("curation_status") == status
                for record in records
            )
            for status in sorted(CURATION_STATUSES)
        }
        if summary.get("disposition_counts") != dispositions:
            errors.append("disposition_counts 不守恒")
        if summary.get("curation_status_counts") != statuses:
            errors.append("curation_status_counts 不守恒")
        if sum(dispositions.values()) != len(records):
            errors.append("record disposition 总数不守恒")
    all_findings = []
    for record in records:
        if isinstance(record, dict):
            all_findings.extend(record.get("findings") or [])
    top_findings = (
        (document.get("errors") or [])
        + (document.get("warnings") or [])
        + (document.get("human_review_required") or [])
    )
    for index, item in enumerate(top_findings):
        errors.extend(validate_finding(item, f"top_findings[{index}]"))
    for key, severity in (
        ("errors", "error"),
        ("warnings", "warning"),
        ("human_review_required", "human_review"),
    ):
        if any(
            isinstance(item, dict) and item.get("severity") != severity
            for item in document.get(key) or []
        ):
            errors.append(f"{key} severity 不一致")
    groups = document.get("duplicate_groups") or []
    known_groups = set()
    for index, group in enumerate(groups):
        if not isinstance(group, dict):
            errors.append(f"duplicate_groups[{index}] 必须是 object")
            continue
        if group.get("view") not in {
            "exact_record",
            "reported_transformation",
            "parent_transformation_candidate",
        }:
            errors.append(f"duplicate_groups[{index}].view 不受控")
        if group.get("automatic_action") != "none":
            errors.append(f"duplicate_groups[{index}] 不得自动操作")
        if len(group.get("record_ids") or []) < 2:
            errors.append(f"duplicate_groups[{index}] 至少两个成员")
        known_groups.add(group.get("group_id"))
    for index, record in enumerate(records):
        if not isinstance(record, dict):
            continue
        unknown = set(record.get("duplicate_memberships") or []) - known_groups
        if unknown:
            errors.append(f"records[{index}] 引用未知 duplicate group")
    for path, key in walk_keys(document):
        if key in FORBIDDEN_KEYS:
            errors.append(f"{path} 是禁止字段")
        if Path(str(document.get(key, ""))).is_absolute() and key.endswith("_path"):
            errors.append(f"{path} 不得保存绝对路径")
    serialized = json.dumps(document, ensure_ascii=False)
    if SECRET_RE.search(serialized):
        errors.append("输出包含疑似凭证")
    lowered = serialized.lower()
    for claim in FORBIDDEN_CLAIMS:
        if claim.lower() in lowered:
            errors.append(f"输出包含禁止结论：{claim}")
    expected_fingerprint = stable_document_fingerprint(document)
    if document.get("result_fingerprint") != expected_fingerprint:
        errors.append("result_fingerprint 不匹配")
    return sorted(set(errors))


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", help="待校验 JSON")
    args = parser.parse_args(argv)
    try:
        document = json.loads(Path(args.input).read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        print(f"ERROR: 无法读取输出：{exc}", file=sys.stderr)
        return 2
    errors = validate_output(document)
    if errors:
        for error in errors:
            print(f"ERROR: {error}", file=sys.stderr)
        return 1
    print("curate-reactions 输出契约校验通过。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
