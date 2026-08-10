#!/usr/bin/env python3
"""Validate search-reactions output and scientific boundaries."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Any, Sequence

from search_reactions import (
    OPERATIONS,
    PROFILE_DEFINITIONS,
    PROVIDERS,
    PROVIDER_STATUSES,
    SCHEMA_VERSION,
    SECRET_RE,
    WORKFLOW,
    sha256_json,
    stable_document_fingerprint,
)

FORBIDDEN_KEYS = {
    "reaction_is_feasible",
    "conditions_are_optimal",
    "safe_to_execute",
    "no_precedent_exists",
    "recommended_conditions",
}
FORBIDDEN_CLAIMS = {
    "反应可行",
    "条件最优",
    "可安全执行",
    "不存在先例",
    "推荐条件",
    "reaction is feasible",
    "conditions are optimal",
    "safe to execute",
    "no precedent exists",
}


def walk(value: Any, path: str = "$") -> list[tuple[str, str, Any]]:
    result = []
    if isinstance(value, dict):
        for key, item in value.items():
            current = f"{path}.{key}"
            result.append((current, str(key), item))
            result.extend(walk(item, current))
    elif isinstance(value, list):
        for index, item in enumerate(value):
            result.extend(walk(item, f"{path}[{index}]"))
    return result


def validate_result(item: Any, path: str, provider: Any) -> list[str]:
    errors = []
    if not isinstance(item, dict):
        return [f"{path} 必须是 object"]
    required = {
        "rank",
        "reaction_id",
        "dataset_id",
        "provider",
        "reaction_smiles",
        "retrieval_mode",
        "fingerprint_profile",
        "raw_score",
        "score_scope",
        "matched_constraints",
        "participants",
        "reported_condition_evidence",
        "yield_measurements",
        "source",
        "license",
        "curation_disposition",
        "quality_findings",
        "result_hash",
    }
    for key in sorted(required):
        if key not in item:
            errors.append(f"{path}.{key} 缺失")
    if item.get("provider") != provider:
        errors.append(f"{path}.provider 与顶层不一致")
    if item.get("curation_disposition") == "rejected":
        errors.append(f"{path} 不得包含 rejected 候选")
    if not isinstance(item.get("rank"), int) or item.get("rank", 0) < 1:
        errors.append(f"{path}.rank 非正整数")
    raw_score = item.get("raw_score")
    if raw_score is not None and (
        not isinstance(raw_score, (int, float))
        or isinstance(raw_score, bool)
        or not 0 <= raw_score <= 1
    ):
        errors.append(f"{path}.raw_score 必须为 null 或 0–1")
    for key in (
        "matched_constraints",
        "participants",
        "yield_measurements",
        "quality_findings",
    ):
        if not isinstance(item.get(key), list):
            errors.append(f"{path}.{key} 必须是 array")
    profile = item.get("fingerprint_profile")
    if item.get("retrieval_mode") == "whole_reaction_similarity":
        if not isinstance(profile, dict):
            errors.append(f"{path}.fingerprint_profile 缺失")
        else:
            profile_id = profile.get("profile_id")
            definition = PROFILE_DEFINITIONS.get(profile_id)
            if definition is None:
                errors.append(f"{path}.fingerprint_profile.profile_id 不受控")
            elif profile.get("metric") != definition["metric"]:
                errors.append(f"{path}.fingerprint_profile.metric 不匹配")
        if item.get("score_scope") != "whole_reaction":
            errors.append(f"{path}.score_scope 必须为 whole_reaction")
    elif profile is not None:
        errors.append(f"{path}.fingerprint_profile 应为 null")
    payload = {
        key: value
        for key, value in item.items()
        if key not in {"rank", "result_hash"}
    }
    if item.get("result_hash") != sha256_json(payload):
        errors.append(f"{path}.result_hash 不匹配")
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
        "operation",
        "provider",
        "provider_status",
        "tool_versions",
        "query_interpretation",
        "options",
        "corpus_summary",
        "results",
        "excluded_records",
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
    if document.get("operation") not in OPERATIONS:
        errors.append("operation 不受控")
    if document.get("provider") not in PROVIDERS:
        errors.append("provider 不受控")
    if document.get("provider_status") not in PROVIDER_STATUSES:
        errors.append("provider_status 不受控")
    versions = document.get("tool_versions")
    if not isinstance(versions, dict):
        errors.append("tool_versions 必须是 object")
    else:
        if versions.get("rdkit") not in {"2025.9.2", "2025.09.2"}:
            errors.append("rdkit 必须固定 2025.9.2")
        if versions.get("ord-schema") != "0.8.3":
            errors.append("ord-schema 必须固定 0.8.3")
    arrays = (
        "results",
        "excluded_records",
        "review_queue",
        "errors",
        "warnings",
        "notices",
    )
    for key in arrays:
        if not isinstance(document.get(key), list):
            errors.append(f"{key} 必须是 array")
    results = document.get("results") or []
    for index, item in enumerate(results):
        errors.extend(
            validate_result(item, f"results[{index}]", document.get("provider"))
        )
    ranks = [item.get("rank") for item in results if isinstance(item, dict)]
    if ranks != list(range(1, len(results) + 1)):
        errors.append("results.rank 必须连续且从 1 开始")
    options = document.get("options")
    if not isinstance(options, dict):
        errors.append("options 必须是 object")
    else:
        top_k = options.get("top_k")
        if not isinstance(top_k, int) or not 1 <= top_k <= 100:
            errors.append("options.top_k 必须为 1–100")
        if len(results) > top_k:
            errors.append("results 数量超过 top_k")
        if document.get("operation") == "search_similar_reactions":
            if (
                document.get("provider_status")
                in {"completed", "completed_zero_hits", "partial"}
                and options.get("fingerprint_profile_id") not in PROFILE_DEFINITIONS
            ):
                errors.append("相似反应检索缺少受控 fingerprint profile")
        elif options.get("fingerprint_profile_id") is not None:
            errors.append("非相似检索不得设置 fingerprint profile")
    status = document.get("provider_status")
    top_errors = document.get("errors") or []
    if status == "completed" and not results:
        errors.append("completed 必须至少有一条结果")
    if status == "completed_zero_hits" and results:
        errors.append("completed_zero_hits 不得有结果")
    if status in {"blocked", "source_timeout", "source_error"} and not top_errors:
        errors.append(f"{status} 必须包含 errors")
    if status in {"completed", "completed_zero_hits"} and top_errors:
        errors.append(f"{status} 不得包含 errors")
    review_ids = {
        item.get("reaction_id")
        for item in document.get("review_queue") or []
        if isinstance(item, dict)
    }
    expected_review = {
        item.get("reaction_id")
        for item in results
        if isinstance(item, dict)
        and item.get("curation_disposition") == "review_required"
    }
    if review_ids != expected_review:
        errors.append("review_queue 与 review_required 结果不一致")
    corpus = document.get("corpus_summary")
    if not isinstance(corpus, dict):
        errors.append("corpus_summary 必须是 object")
    else:
        for key in ("input_records", "searchable_records", "excluded_records"):
            if not isinstance(corpus.get(key), int) or corpus.get(key, -1) < 0:
                errors.append(f"corpus_summary.{key} 必须是非负整数")
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
