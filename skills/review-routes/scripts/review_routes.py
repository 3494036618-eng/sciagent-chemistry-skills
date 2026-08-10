#!/usr/bin/env python3
"""Review existing synthesis routes with deterministic evidence contracts."""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import re
import sys
import time
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Sequence

SCHEMA_VERSION = "1.0.0"
WORKFLOW = "review-routes"
RULESET_VERSION = "1.0.0"
INPUT_PROFILES = {
    "normalized_route_v1",
    "aizynthfinder_json",
    "paroutes_v2_json",
}
REVIEW_STATUSES = {"completed", "partial", "not_run", "error"}
DISPOSITIONS = {"ready_for_expert_review", "review_required", "blocked"}
PRECEDENT_LEVELS = {
    "exact_record",
    "exact_transformation",
    "similar_reaction",
    "component_only",
    "completed_zero_hits",
    "source_timeout",
    "source_error",
    "not_run",
}
MAX_ROUTES = 20
MAX_STEPS_PER_ROUTE = 50
MAX_TOTAL_NODES = 5000
MAX_STEP_ARTIFACTS = MAX_ROUTES * MAX_STEPS_PER_ROUTE
TEMPORAL_KEYS = {
    "generated_at_utc",
    "runtime_seconds",
    "retrieved_at_utc",
    "elapsed_seconds",
    "result_fingerprint",
}
SECRET_RE = re.compile(
    r"(?-i:ark-)[A-Za-z0-9_-]{12,}|"
    r"Bearer\s+[A-Za-z0-9._~+/=-]{12,}|"
    r"(?:Authorization|Cookie|Token|Api[_ -]?Key)\s*[:=]\s*\S{12,}",
    re.IGNORECASE,
)

RULE_MESSAGES = {
    "E-INPUT-SCHEMA-001": "输入字段、枚举或显式选项不符合冻结合同。",
    "E-INPUT-HASH-001": "来源或路线缺少有效 SHA-256。",
    "E-RESOURCE-LIMIT-001": "路线、步骤、节点或 artifact 超过首版资源上限。",
    "E-ROUTE-ID-001": "route_id 缺失或批内重复。",
    "E-ROUTE-TOPOLOGY-001": "路线不是单根、连通、无环的 mol/reaction 交替树。",
    "E-MOLECULE-STRUCTURE-001": "路线分子结构缺失或 RDKit 无法解析。",
    "E-STEP-REACTION-001": "路线步骤无法形成与父产品一致的单步反应。",
    "E-ARTIFACT-FINGERPRINT-001": "上游 artifact 指纹缺失、类型错误或内容不匹配。",
    "E-STEP-HASH-MISMATCH-001": "step artifact 未绑定到当前 step reaction hash。",
    "E-CURATION-REJECTED-001": "步骤被 curate-reactions 拒绝，路线不得进入 ready 状态。",
    "E-PROFILE-MISMATCH-001": "同一步先例 artifact 混入多个不可比 fingerprint profile。",
    "E-PICKLE-INPUT-001": "禁止读取 pickle 或其他可执行反序列化格式。",
    "W-TARGET-MISMATCH-001": "请求目标与路线根分子不一致。",
    "W-CURATION-REVIEW-001": "步骤继承 curate-reactions 人工复核状态。",
    "W-PRECEDENT-SIMILAR-001": "步骤只有相似反应先例，不证明该步骤可行。",
    "W-PRECEDENT-COMPONENT-001": "步骤只有组分级命中，不等于转化先例。",
    "W-PRECEDENT-ZERO-001": "当前 provider/query 返回 0 hit，不代表不存在先例。",
    "W-PRECEDENT-TIMEOUT-001": "先例 provider 超时，不能解释为 0 hit。",
    "W-PRECEDENT-ERROR-001": "先例 provider 错误，不能解释为 0 hit。",
    "W-PRECEDENT-NOT-RUN-001": "步骤尚未执行先例检索。",
    "W-INVENTORY-MISSING-001": "缺少可审计库存快照或终端前体状态。",
    "W-INVENTORY-LICENSE-001": "库存快照缺少明确许可。",
    "W-SOURCE-LICENSE-001": "路线或先例来源缺少明确许可。",
    "W-CONSTRAINT-VIOLATION-001": "路线违反用户显式提供的项目约束。",
    "W-ROUTE-DUPLICATE-001": "路线与其他候选具有相同结构签名，仅分组不删除。",
}


class DependencyFailure(RuntimeError):
    """Fixed chemistry dependency is unavailable."""


class InputFailure(ValueError):
    """Input cannot be processed under the frozen contract."""


class ResourceFailure(InputFailure):
    """Input exceeds a frozen resource limit."""


def now_utc() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def canonical_json(value: Any) -> str:
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"))


def sha256_json(value: Any) -> str:
    return hashlib.sha256(canonical_json(value).encode("utf-8")).hexdigest()


def valid_sha256(value: Any) -> bool:
    return isinstance(value, str) and re.fullmatch(r"[0-9a-f]{64}", value) is not None


def without_temporal(value: Any) -> Any:
    if isinstance(value, dict):
        return {
            key: without_temporal(item)
            for key, item in value.items()
            if key not in TEMPORAL_KEYS
        }
    if isinstance(value, list):
        return [without_temporal(item) for item in value]
    return value


def stable_document_fingerprint(document: dict[str, Any]) -> str:
    return sha256_json(without_temporal(document))


def artifact_fingerprint(artifact: dict[str, Any]) -> str:
    workflow = artifact.get("workflow")
    if workflow == "curate-reactions":
        return sha256_json(
            {
                key: value
                for key, value in artifact.items()
                if key
                not in {"generated_at_utc", "runtime_seconds", "result_fingerprint"}
            }
        )
    return stable_document_fingerprint(artifact)


def load_toolkit() -> dict[str, Any]:
    try:
        import rdkit
        from rdkit import Chem, rdBase
    except ImportError as error:
        raise DependencyFailure(
            "需要 rdkit==2025.9.2；请在隔离环境安装 scripts/requirements.txt。"
        ) from error
    if rdkit.__version__ not in {"2025.9.2", "2025.09.2"}:
        raise DependencyFailure(f"需要 rdkit==2025.9.2，当前为 {rdkit.__version__}。")
    return {"rdkit": rdkit, "Chem": Chem, "rdBase": rdBase}


def tool_versions(toolkit: dict[str, Any]) -> dict[str, str]:
    return {
        "python": platform.python_version(),
        "rdkit": toolkit["rdkit"].__version__,
        "review-routes": RULESET_VERSION,
    }


def finding(
    code: str,
    severity: str,
    field_path: str,
    *,
    evidence: Any = None,
) -> dict[str, Any]:
    item = {
        "code": code,
        "severity": severity,
        "field_path": field_path,
        "message": RULE_MESSAGES[code],
        "evidence": [] if evidence is None else evidence if isinstance(evidence, list) else [evidence],
    }
    return item


def parse_molecule(value: Any, toolkit: dict[str, Any]) -> tuple[Any | None, str | None]:
    if not isinstance(value, str) or not value.strip():
        return None, None
    try:
        with toolkit["rdBase"].BlockLogs():
            molecule = toolkit["Chem"].MolFromSmiles(value)
    except Exception:
        molecule = None
    if molecule is None:
        return None, None
    for atom in molecule.GetAtoms():
        atom.SetAtomMapNum(0)
    canonical = toolkit["Chem"].MolToSmiles(
        molecule, canonical=True, isomericSmiles=True
    )
    return molecule, canonical


def split_reaction_smiles(value: Any) -> tuple[list[str], list[str], list[str]]:
    if not isinstance(value, str) or not value.strip():
        raise InputFailure("reaction SMILES 为空。")
    text = value.strip()
    if ">>" in text:
        if text.count(">>") != 1:
            raise InputFailure("reaction SMILES 不是单步两段或三段形式。")
        left, right = text.split(">>")
        middle = ""
    else:
        parts = text.split(">")
        if len(parts) != 3:
            raise InputFailure("reaction SMILES 不是单步两段或三段形式。")
        left, middle, right = parts
    inputs = [item for item in left.split(".") if item]
    agents = [item for item in middle.split(".") if item]
    outputs = [item for item in right.split(".") if item]
    if not inputs or not outputs:
        raise InputFailure("reaction SMILES 缺少输入或输出。")
    return inputs, agents, outputs


def canonical_reaction_smiles(value: str, toolkit: dict[str, Any]) -> str:
    inputs, agents, outputs = split_reaction_smiles(value)

    def canonical_side(side: list[str]) -> list[str]:
        result = []
        for structure in side:
            _, canonical = parse_molecule(structure, toolkit)
            if canonical is None:
                raise InputFailure(f"reaction component 无法解析：{structure!r}")
            result.append(canonical)
        return sorted(result)

    return ">".join(
        (
            ".".join(canonical_side(inputs)),
            ".".join(canonical_side(agents)),
            ".".join(canonical_side(outputs)),
        )
    )


def node_kind(node: Any) -> str | None:
    if not isinstance(node, dict):
        return None
    value = str(node.get("type") or node.get("kind") or "").strip().lower()
    if value in {"mol", "molecule"}:
        return "mol"
    if value in {"reaction", "rxn"} or node.get("is_reaction") is True:
        return "reaction"
    return None


def child_nodes(node: dict[str, Any]) -> list[Any]:
    value = node.get("children")
    return value if isinstance(value, list) else []


def reaction_text(node: dict[str, Any]) -> str | None:
    metadata = node.get("metadata")
    if not isinstance(metadata, dict):
        metadata = {}
    for value in (
        node.get("reaction_smiles"),
        metadata.get("reaction_smiles"),
        metadata.get("mapped_reaction"),
        metadata.get("smiles"),
        metadata.get("rsmi"),
        node.get("smiles"),
    ):
        if isinstance(value, str) and (">>" in value or value.count(">") == 2):
            return value
    return None


def normalize_routes(
    request: dict[str, Any],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    profile = request.get("input_profile")
    values = request.get("routes")
    if not isinstance(values, list):
        raise InputFailure("routes 必须是 array。")
    if len(values) > MAX_ROUTES:
        raise ResourceFailure(f"routes 不得超过 {MAX_ROUTES} 条。")
    normalized = []
    errors = []
    seen: set[str] = set()
    for index, value in enumerate(values):
        if not isinstance(value, dict):
            errors.append(
                finding(
                    "E-INPUT-SCHEMA-001",
                    "error",
                    f"routes[{index}]",
                )
            )
            continue
        if profile == "paroutes_v2_json":
            tree = value
            route_id = value.get("route_id") or f"paroutes-{sha256_json(value)[:16]}"
            backend = "paroutes"
            backend_rank = index + 1
            backend_score = None
        else:
            tree = value.get("tree")
            route_id = value.get("route_id")
            if not route_id and profile == "aizynthfinder_json":
                route_id = f"aizynth-{sha256_json(value)[:16]}"
            backend = value.get("backend") or (
                "aizynthfinder" if profile == "aizynthfinder_json" else None
            )
            backend_rank = value.get("backend_rank", value.get("rank"))
            backend_score = value.get("backend_score", value.get("score"))
        if not isinstance(route_id, str) or not route_id or route_id in seen:
            errors.append(
                finding(
                    "E-ROUTE-ID-001",
                    "error",
                    f"routes[{index}].route_id",
                    evidence=route_id,
                )
            )
            continue
        seen.add(route_id)
        normalized.append(
            {
                "route_id": route_id,
                "backend": backend,
                "backend_rank": backend_rank,
                "backend_score": backend_score,
                "source_route_hash": sha256_json(value),
                "tree": tree,
                "source_index": index,
            }
        )
    return normalized, errors


def analyze_route_tree(
    route: dict[str, Any],
    toolkit: dict[str, Any],
) -> dict[str, Any]:
    route_id = route["route_id"]
    findings: list[dict[str, Any]] = []
    steps: list[dict[str, Any]] = []
    leaves: list[dict[str, Any]] = []
    node_count = 0
    active: set[int] = set()

    def walk_molecule(node: Any, path: tuple[int, ...]) -> tuple[str | None, int]:
        nonlocal node_count
        node_count += 1
        field_path = f"routes[{route['source_index']}].tree" + "".join(
            f".children[{index}]" for index in path
        )
        if not isinstance(node, dict) or node_kind(node) != "mol":
            findings.append(finding("E-ROUTE-TOPOLOGY-001", "error", field_path))
            return None, 0
        identity = id(node)
        if identity in active:
            findings.append(
                finding(
                    "E-ROUTE-TOPOLOGY-001",
                    "error",
                    field_path,
                    evidence="cycle",
                )
            )
            return None, 0
        active.add(identity)
        _, canonical = parse_molecule(node.get("smiles"), toolkit)
        if canonical is None:
            findings.append(
                finding(
                    "E-MOLECULE-STRUCTURE-001",
                    "error",
                    f"{field_path}.smiles",
                    evidence=node.get("smiles"),
                )
            )
        children = child_nodes(node)
        if not children:
            leaves.append(
                {
                    "structure": canonical,
                    "reported_structure": node.get("smiles"),
                    "reported_in_stock": (
                        node.get("in_stock")
                        if isinstance(node.get("in_stock"), bool)
                        else None
                    ),
                    "path": list(path),
                }
            )
            active.remove(identity)
            return canonical, 0
        if len(children) != 1 or node_kind(children[0]) != "reaction":
            findings.append(
                finding(
                    "E-ROUTE-TOPOLOGY-001",
                    "error",
                    f"{field_path}.children",
                    evidence=f"expected one reaction child, got {len(children)}",
                )
            )
            active.remove(identity)
            return canonical, 0
        reaction = children[0]
        node_count += 1
        reaction_path = (*path, 0)
        reaction_field = f"{field_path}.children[0]"
        reaction_identity = id(reaction)
        if reaction_identity in active:
            findings.append(
                finding(
                    "E-ROUTE-TOPOLOGY-001",
                    "error",
                    reaction_field,
                    evidence="cycle",
                )
            )
            active.remove(identity)
            return canonical, 0
        active.add(reaction_identity)
        precursors = child_nodes(reaction)
        if not precursors or any(node_kind(item) != "mol" for item in precursors):
            findings.append(
                finding(
                    "E-ROUTE-TOPOLOGY-001",
                    "error",
                    f"{reaction_field}.children",
                    evidence="reaction requires molecule precursor children",
                )
            )
        precursor_values = []
        child_depths = []
        for child_index, precursor in enumerate(precursors):
            value, depth = walk_molecule(precursor, (*reaction_path, child_index))
            if value:
                precursor_values.append(value)
            child_depths.append(depth)
        reported_reaction = reaction_text(reaction)
        if reported_reaction is None and canonical and precursor_values:
            reported_reaction = f"{'.'.join(precursor_values)}>>{canonical}"
        canonical_reaction = None
        agents: list[str] = []
        if reported_reaction:
            try:
                canonical_reaction = canonical_reaction_smiles(
                    reported_reaction, toolkit
                )
                inputs, agents, outputs = split_reaction_smiles(canonical_reaction)
                if canonical not in outputs:
                    findings.append(
                        finding(
                            "E-STEP-REACTION-001",
                            "error",
                            reaction_field,
                            evidence="reaction output does not match parent product",
                        )
                    )
                available = set(inputs) | set(agents)
                missing = sorted(set(precursor_values) - available)
                if missing:
                    findings.append(
                        finding(
                            "E-STEP-REACTION-001",
                            "error",
                            reaction_field,
                            evidence={"precursors_missing_from_reaction": missing},
                        )
                    )
            except InputFailure as error:
                findings.append(
                    finding(
                        "E-STEP-REACTION-001",
                        "error",
                        reaction_field,
                        evidence=str(error),
                    )
                )
        else:
            findings.append(
                finding("E-STEP-REACTION-001", "error", reaction_field)
            )
        step_hash = (
            hashlib.sha256(canonical_reaction.encode("utf-8")).hexdigest()
            if canonical_reaction
            else None
        )
        step_id = "step-" + hashlib.sha256(
            f"{route_id}:{reaction_path}:{step_hash}".encode("utf-8")
        ).hexdigest()[:16]
        metadata = reaction.get("metadata")
        steps.append(
            {
                "step_id": step_id,
                "step_reaction_hash": step_hash,
                "path": list(reaction_path),
                "reported_reaction": reported_reaction,
                "canonical_reaction": canonical_reaction,
                "product": canonical,
                "precursors": sorted(precursor_values),
                "agents": agents,
                "backend_metadata": metadata if isinstance(metadata, dict) else {},
            }
        )
        active.remove(reaction_identity)
        active.remove(identity)
        return canonical, 1 + max(child_depths, default=0)

    tree = route.get("tree")
    if node_kind(tree) != "mol":
        findings.append(
            finding(
                "E-ROUTE-TOPOLOGY-001",
                "error",
                f"routes[{route['source_index']}].tree",
                evidence="root must be a molecule",
            )
        )
        root_structure, longest = None, 0
    else:
        root_structure, longest = walk_molecule(tree, ())
    if len(steps) > MAX_STEPS_PER_ROUTE:
        findings.append(
            finding(
                "E-RESOURCE-LIMIT-001",
                "error",
                f"routes[{route['source_index']}].tree",
                evidence={"steps": len(steps)},
            )
        )
    return {
        "root_structure": root_structure,
        "steps": steps,
        "leaves": leaves,
        "node_count": node_count,
        "longest_linear_sequence": longest,
        "branch_count": sum(max(0, len(step["precursors"]) - 1) for step in steps),
        "findings": findings,
    }


def validate_artifact(artifact: Any, workflow: str) -> tuple[dict[str, Any] | None, str | None]:
    if not isinstance(artifact, dict) or artifact.get("workflow") != workflow:
        return None, "workflow_or_type"
    actual = artifact.get("result_fingerprint")
    if not valid_sha256(actual) or actual != artifact_fingerprint(artifact):
        return None, "fingerprint"
    return artifact, None


def curation_evidence(
    artifact: Any,
    step: dict[str, Any],
    toolkit: dict[str, Any],
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    result = {
        "status": "not_run",
        "disposition": None,
        "findings": [],
        "artifact_fingerprint": None,
    }
    findings = []
    if artifact is None:
        return result, findings
    value, error = validate_artifact(artifact, "curate-reactions")
    if error:
        findings.append(
            finding(
                "E-ARTIFACT-FINGERPRINT-001",
                "error",
                "step_artifacts.curation_artifact",
                evidence=error,
            )
        )
        return result, findings
    records = value.get("records")
    if not isinstance(records, list):
        findings.append(
            finding(
                "E-ARTIFACT-FINGERPRINT-001",
                "error",
                "step_artifacts.curation_artifact.records",
            )
        )
        return result, findings
    match = None
    for record in records:
        if not isinstance(record, dict):
            continue
        reaction_value = record.get("reaction_smiles")
        if isinstance(reaction_value, dict):
            reaction_value = reaction_value.get("canonical_unmapped") or reaction_value.get(
                "reported"
            )
        if not isinstance(reaction_value, str):
            continue
        try:
            canonical = canonical_reaction_smiles(reaction_value, toolkit)
        except InputFailure:
            continue
        if (
            hashlib.sha256(canonical.encode("utf-8")).hexdigest()
            == step["step_reaction_hash"]
        ):
            match = record
            break
    if match is None:
        findings.append(
            finding(
                "E-STEP-HASH-MISMATCH-001",
                "error",
                "step_artifacts.curation_artifact",
                evidence="no record matches step reaction hash",
            )
        )
        return result, findings
    disposition = match.get("disposition")
    status = match.get("curation_status")
    result = {
        "status": status,
        "disposition": disposition,
        "findings": match.get("findings") if isinstance(match.get("findings"), list) else [],
        "artifact_fingerprint": value.get("result_fingerprint"),
    }
    if disposition == "rejected":
        findings.append(
            finding(
                "E-CURATION-REJECTED-001",
                "error",
                "step_artifacts.curation_artifact",
            )
        )
    elif disposition == "review_required":
        findings.append(
            finding(
                "W-CURATION-REVIEW-001",
                "warning",
                "step_artifacts.curation_artifact",
            )
        )
    return result, findings


def precedent_level(document: dict[str, Any]) -> str:
    status = document.get("provider_status")
    if status in {"completed_zero_hits", "source_timeout", "source_error"}:
        return status
    results = document.get("results")
    if not isinstance(results, list) or not results:
        return "not_run"
    levels = []
    for result in results:
        if not isinstance(result, dict):
            continue
        mode = result.get("retrieval_mode")
        if mode == "exact_id":
            levels.append("exact_record")
        elif mode == "reaction_smarts":
            levels.append("exact_transformation")
        elif mode == "whole_reaction_similarity":
            levels.append("similar_reaction")
        elif mode == "component_and_filter":
            levels.append("component_only")
    order = [
        "exact_record",
        "exact_transformation",
        "similar_reaction",
        "component_only",
    ]
    return next((item for item in order if item in levels), "not_run")


def precedent_evidence(
    artifact: Any,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    result = {
        "provider_status": "not_run",
        "match_level": "not_run",
        "profile_ids": [],
        "reported_condition_evidence": [],
        "reported_yield_evidence": [],
        "sources": [],
        "licenses": [],
        "artifact_fingerprint": None,
    }
    findings = []
    if artifact is None:
        findings.append(
            finding(
                "W-PRECEDENT-NOT-RUN-001",
                "warning",
                "step_artifacts.precedent_artifact",
            )
        )
        return result, findings
    value, error = validate_artifact(artifact, "search-reactions")
    if error:
        findings.append(
            finding(
                "E-ARTIFACT-FINGERPRINT-001",
                "error",
                "step_artifacts.precedent_artifact",
                evidence=error,
            )
        )
        return result, findings
    level = precedent_level(value)
    results = value.get("results") if isinstance(value.get("results"), list) else []
    profiles = sorted(
        {
            profile.get("profile_id")
            for item in results
            if isinstance(item, dict)
            and isinstance((profile := item.get("fingerprint_profile")), dict)
            and profile.get("profile_id")
        }
    )
    conditions = [
        item.get("reported_condition_evidence")
        for item in results
        if isinstance(item, dict) and item.get("reported_condition_evidence")
    ]
    yields = [
        item.get("yield_measurements")
        for item in results
        if isinstance(item, dict) and item.get("yield_measurements")
    ]
    sources = [
        item.get("source")
        for item in results
        if isinstance(item, dict) and isinstance(item.get("source"), dict)
    ]
    licenses = sorted(
        {
            str(item.get("license"))
            for item in results
            if isinstance(item, dict) and item.get("license")
        }
    )
    result = {
        "provider_status": value.get("provider_status"),
        "match_level": level,
        "profile_ids": profiles,
        "reported_condition_evidence": conditions,
        "reported_yield_evidence": yields,
        "sources": sources,
        "licenses": licenses,
        "artifact_fingerprint": value.get("result_fingerprint"),
    }
    code_by_level = {
        "similar_reaction": "W-PRECEDENT-SIMILAR-001",
        "component_only": "W-PRECEDENT-COMPONENT-001",
        "completed_zero_hits": "W-PRECEDENT-ZERO-001",
        "source_timeout": "W-PRECEDENT-TIMEOUT-001",
        "source_error": "W-PRECEDENT-ERROR-001",
        "not_run": "W-PRECEDENT-NOT-RUN-001",
    }
    if level in code_by_level:
        findings.append(
            finding(
                code_by_level[level],
                "warning",
                "step_artifacts.precedent_artifact",
            )
        )
    if not licenses:
        findings.append(
            finding(
                "W-SOURCE-LICENSE-001",
                "warning",
                "step_artifacts.precedent_artifact.results.license",
            )
        )
    if len(profiles) > 1:
        findings.append(
            finding(
                "E-PROFILE-MISMATCH-001",
                "error",
                "step_artifacts.precedent_artifact.results.fingerprint_profile",
                evidence=profiles,
            )
        )
    return result, findings


def build_artifact_index(
    entries: Any,
) -> tuple[dict[tuple[str, str], dict[str, Any]], list[dict[str, Any]]]:
    if entries is None:
        return {}, []
    if not isinstance(entries, list):
        return {}, [finding("E-INPUT-SCHEMA-001", "error", "step_artifacts")]
    if len(entries) > MAX_STEP_ARTIFACTS:
        return {}, [
            finding("E-RESOURCE-LIMIT-001", "error", "step_artifacts")
        ]
    index = {}
    errors = []
    for position, entry in enumerate(entries):
        if not isinstance(entry, dict):
            errors.append(
                finding(
                    "E-INPUT-SCHEMA-001",
                    "error",
                    f"step_artifacts[{position}]",
                )
            )
            continue
        key = (entry.get("route_id"), entry.get("step_id"))
        if not all(isinstance(item, str) and item for item in key) or key in index:
            errors.append(
                finding(
                    "E-INPUT-SCHEMA-001",
                    "error",
                    f"step_artifacts[{position}]",
                    evidence="invalid or duplicate route_id/step_id",
                )
            )
            continue
        index[key] = entry
    return index, errors


def normalize_inventory(
    value: Any,
    toolkit: dict[str, Any],
) -> tuple[dict[str, Any] | None, dict[str, dict[str, Any]], list[dict[str, Any]]]:
    if value is None:
        return None, {}, [
            finding("W-INVENTORY-MISSING-001", "warning", "inventory_snapshot")
        ]
    if not isinstance(value, dict):
        return None, {}, [
            finding("E-INPUT-SCHEMA-001", "error", "inventory_snapshot")
        ]
    required = ("snapshot_id", "captured_at_utc", "source", "records")
    if any(not value.get(key) for key in required) or not isinstance(
        value.get("records"), list
    ):
        return value, {}, [
            finding("W-INVENTORY-MISSING-001", "warning", "inventory_snapshot")
        ]
    findings = []
    if not value.get("license"):
        findings.append(
            finding(
                "W-INVENTORY-LICENSE-001",
                "warning",
                "inventory_snapshot.license",
            )
        )
    index = {}
    for position, record in enumerate(value["records"]):
        if not isinstance(record, dict):
            continue
        _, canonical = parse_molecule(record.get("structure"), toolkit)
        status = record.get("status")
        if canonical and status in {"in_stock", "not_in_stock", "unknown"}:
            index[canonical] = {
                "status": status,
                "source_record": record,
                "position": position,
            }
    return value, index, findings


def normalize_constraints(value: Any, toolkit: dict[str, Any]) -> dict[str, Any]:
    if value is None:
        return {}
    if not isinstance(value, dict):
        raise InputFailure("constraints 必须是 object。")
    allowed = {
        "max_steps",
        "max_precursors",
        "require_all_leaves_in_stock",
        "minimum_exact_or_transformation_coverage",
        "forbidden_starting_materials",
    }
    unknown = sorted(set(value) - allowed)
    if unknown:
        raise InputFailure(f"不支持的 constraints：{unknown}")
    result = {}
    for key in ("max_steps", "max_precursors"):
        if key in value:
            item = value[key]
            if not isinstance(item, int) or isinstance(item, bool) or item < 0:
                raise InputFailure(f"constraints.{key} 必须是非负整数。")
            result[key] = item
    key = "minimum_exact_or_transformation_coverage"
    if key in value:
        item = value[key]
        if (
            not isinstance(item, (int, float))
            or isinstance(item, bool)
            or not 0 <= item <= 1
        ):
            raise InputFailure(f"constraints.{key} 必须为 0–1。")
        result[key] = float(item)
    if "require_all_leaves_in_stock" in value:
        if not isinstance(value["require_all_leaves_in_stock"], bool):
            raise InputFailure("require_all_leaves_in_stock 必须是 boolean。")
        result["require_all_leaves_in_stock"] = value[
            "require_all_leaves_in_stock"
        ]
    if "forbidden_starting_materials" in value:
        values = value["forbidden_starting_materials"]
        if not isinstance(values, list) or not all(
            isinstance(item, str) for item in values
        ):
            raise InputFailure("forbidden_starting_materials 必须是 string array。")
        canonical_values = []
        for item in values:
            _, canonical = parse_molecule(item, toolkit)
            if canonical is None:
                raise InputFailure(f"禁用前体无法解析：{item!r}")
            canonical_values.append(canonical)
        result["forbidden_starting_materials"] = sorted(set(canonical_values))
    return result


def route_signature(route: dict[str, Any]) -> str:
    payload = {
        "target": route.get("target_structure"),
        "steps": sorted(
            (
                {
                    "reaction": item.get("canonical_reaction"),
                    "product": item.get("product"),
                    "precursors": item.get("precursors"),
                }
                for item in route.get("step_reviews") or []
            ),
            key=canonical_json,
        ),
        "leaves": sorted(
            item.get("structure")
            for item in route.get("terminal_precursors") or []
            if item.get("structure")
        ),
    }
    return "route:" + sha256_json(payload)[:24]


def finalize_route(route: dict[str, Any]) -> None:
    findings = route["findings"]
    severities = {item["severity"] for item in findings}
    codes = sorted({item["code"] for item in findings})
    if "error" in severities:
        route["review_status"] = "error"
        route["disposition"] = "blocked"
    else:
        partial_codes = {
            "W-CURATION-REVIEW-001",
            "W-PRECEDENT-TIMEOUT-001",
            "W-PRECEDENT-ERROR-001",
            "W-PRECEDENT-NOT-RUN-001",
        }
        route["review_status"] = (
            "partial" if any(code in partial_codes for code in codes) else "completed"
        )
        route["disposition"] = (
            "review_required" if findings else "ready_for_expert_review"
        )
    route["human_review_required"] = [
        code for code in codes if code.startswith("W-")
    ]


def process_route(
    route: dict[str, Any],
    *,
    target_structure: str | None,
    artifact_index: dict[tuple[str, str], dict[str, Any]],
    inventory: dict[str, Any] | None,
    inventory_index: dict[str, dict[str, Any]],
    inventory_findings: list[dict[str, Any]],
    constraints: dict[str, Any],
    source_license: Any,
    toolkit: dict[str, Any],
) -> dict[str, Any]:
    analysis = analyze_route_tree(route, toolkit)
    findings = list(analysis["findings"])
    if target_structure and analysis["root_structure"] != target_structure:
        findings.append(
            finding(
                "W-TARGET-MISMATCH-001",
                "warning",
                f"routes[{route['source_index']}].tree.smiles",
                evidence={
                    "request": target_structure,
                    "route": analysis["root_structure"],
                },
            )
        )
    if not source_license:
        findings.append(
            finding("W-SOURCE-LICENSE-001", "warning", "source.license")
        )
    step_reviews = []
    coverage = Counter()
    for step in analysis["steps"]:
        entry = artifact_index.get((route["route_id"], step["step_id"]))
        step_findings = []
        if entry is not None and entry.get("step_reaction_hash") != step[
            "step_reaction_hash"
        ]:
            step_findings.append(
                finding(
                    "E-STEP-HASH-MISMATCH-001",
                    "error",
                    "step_artifacts.step_reaction_hash",
                    evidence={
                        "expected": step["step_reaction_hash"],
                        "actual": entry.get("step_reaction_hash"),
                    },
                )
            )
        curation, curation_findings = curation_evidence(
            entry.get("curation_artifact") if entry else None,
            step,
            toolkit,
        )
        precedent, precedent_findings = precedent_evidence(
            entry.get("precedent_artifact") if entry else None
        )
        step_findings.extend(curation_findings)
        step_findings.extend(precedent_findings)
        coverage[precedent["match_level"]] += 1
        findings.extend(
            {
                **item,
                "evidence": [
                    {"step_id": step["step_id"]},
                    *item.get("evidence", []),
                ],
            }
            for item in step_findings
        )
        step_reviews.append(
            {
                **step,
                "curation": curation,
                "precedent": precedent,
                "findings": step_findings,
                "review_required": sorted(
                    {item["code"] for item in step_findings}
                ),
            }
        )
    terminal_precursors = []
    for leaf in analysis["leaves"]:
        inventory_record = inventory_index.get(leaf["structure"])
        if inventory_record:
            status = inventory_record["status"]
            source = "inventory_snapshot"
        elif leaf["reported_in_stock"] is True:
            status = "reported_in_stock"
            source = "route_export"
        elif leaf["reported_in_stock"] is False:
            status = "reported_not_in_stock"
            source = "route_export"
        else:
            status = "unknown"
            source = None
        terminal_precursors.append(
            {
                **leaf,
                "inventory_status": status,
                "inventory_source": source,
            }
        )
    findings.extend(inventory_findings)
    if not inventory_index or any(
        item["inventory_status"] in {"unknown", "reported_not_in_stock"}
        for item in terminal_precursors
    ):
        findings.append(
            finding(
                "W-INVENTORY-MISSING-001",
                "warning",
                "inventory_snapshot",
            )
        )
    step_count = len(step_reviews)
    exact_count = coverage["exact_record"] + coverage["exact_transformation"]
    exact_coverage = exact_count / step_count if step_count else 0.0
    constraint_results = []

    def constraint(name: str, passed: bool, observed: Any, expected: Any) -> None:
        constraint_results.append(
            {
                "constraint": name,
                "passed": passed,
                "observed": observed,
                "expected": expected,
            }
        )
        if not passed:
            findings.append(
                finding(
                    "W-CONSTRAINT-VIOLATION-001",
                    "warning",
                    f"constraints.{name}",
                    evidence={"observed": observed, "expected": expected},
                )
            )

    if "max_steps" in constraints:
        constraint("max_steps", step_count <= constraints["max_steps"], step_count, constraints["max_steps"])
    if "max_precursors" in constraints:
        constraint(
            "max_precursors",
            len(terminal_precursors) <= constraints["max_precursors"],
            len(terminal_precursors),
            constraints["max_precursors"],
        )
    if constraints.get("require_all_leaves_in_stock"):
        all_in_stock = bool(terminal_precursors) and all(
            item["inventory_status"] in {"in_stock", "reported_in_stock"}
            for item in terminal_precursors
        )
        constraint("require_all_leaves_in_stock", all_in_stock, all_in_stock, True)
    if "minimum_exact_or_transformation_coverage" in constraints:
        constraint(
            "minimum_exact_or_transformation_coverage",
            exact_coverage
            >= constraints["minimum_exact_or_transformation_coverage"],
            exact_coverage,
            constraints["minimum_exact_or_transformation_coverage"],
        )
    forbidden = set(constraints.get("forbidden_starting_materials") or [])
    if forbidden:
        present = sorted(
            {
                item["structure"]
                for item in terminal_precursors
                if item["structure"] in forbidden
            }
        )
        constraint("forbidden_starting_materials", not present, present, [])
    output = {
        "route_id": route["route_id"],
        "source_route_hash": route["source_route_hash"],
        "backend_metadata": {
            "backend": route.get("backend"),
            "backend_rank": route.get("backend_rank"),
            "backend_score": route.get("backend_score"),
        },
        "route_signature": None,
        "target_structure": analysis["root_structure"],
        "topology_status": (
            "invalid"
            if any(item["severity"] == "error" for item in analysis["findings"])
            else "valid"
        ),
        "node_count": analysis["node_count"],
        "step_count": step_count,
        "longest_linear_sequence": analysis["longest_linear_sequence"],
        "branch_count": analysis["branch_count"],
        "terminal_precursors": terminal_precursors,
        "inventory_snapshot": (
            {
                key: inventory.get(key)
                for key in ("snapshot_id", "captured_at_utc", "source", "license")
            }
            if isinstance(inventory, dict)
            else None
        ),
        "inventory_coverage": (
            sum(
                item["inventory_status"] in {"in_stock", "not_in_stock"}
                for item in terminal_precursors
            )
            / len(terminal_precursors)
            if terminal_precursors
            else 0.0
        ),
        "precedent_coverage_by_level": {
            level: coverage.get(level, 0) for level in sorted(PRECEDENT_LEVELS)
        },
        "exact_or_transformation_coverage": exact_coverage,
        "weakest_steps": [
            item["step_id"]
            for item in step_reviews
            if item["precedent"]["match_level"]
            not in {"exact_record", "exact_transformation"}
            or item["curation"]["disposition"] != "ready_for_search"
        ],
        "constraint_results": constraint_results,
        "step_reviews": step_reviews,
        "findings": findings,
        "review_status": None,
        "disposition": None,
        "human_review_required": [],
        "duplicate_memberships": [],
    }
    output["route_signature"] = route_signature(output)
    finalize_route(output)
    return output


def duplicate_groups(routes: list[dict[str, Any]]) -> list[dict[str, Any]]:
    groups = defaultdict(list)
    for route in routes:
        groups[route["route_signature"]].append(route)
    output = []
    for signature, members in sorted(groups.items()):
        if len(members) < 2:
            continue
        group_id = "duplicate-" + hashlib.sha256(signature.encode("utf-8")).hexdigest()[:12]
        route_ids = sorted(item["route_id"] for item in members)
        output.append(
            {
                "group_id": group_id,
                "route_signature": signature,
                "route_ids": route_ids,
                "member_count": len(route_ids),
            }
        )
        for route in members:
            route["duplicate_memberships"].append(group_id)
            route["findings"].append(
                finding(
                    "W-ROUTE-DUPLICATE-001",
                    "warning",
                    "route_signature",
                    evidence={"group_id": group_id, "route_ids": route_ids},
                )
            )
            finalize_route(route)
    return output


def process_request(
    request: dict[str, Any],
    *,
    generated_at_utc: str | None = None,
) -> dict[str, Any]:
    started = time.perf_counter()
    toolkit = load_toolkit()
    top_findings: list[dict[str, Any]] = []
    source = request.get("source")
    source = source if isinstance(source, dict) else {}
    public_source = {
        key: source.get(key)
        for key in ("identifier", "content_sha256", "license")
    }
    options = request.get("options")
    if not isinstance(options, dict):
        options = {}
    normalized_options = {
        "comparison_mode": options.get("comparison_mode"),
        "preserve_backend_order": options.get("preserve_backend_order"),
        "automatic_route_ranking": False,
        "network_access": False,
        "pickle_allowed": False,
    }
    routes: list[dict[str, Any]] = []
    routes_fingerprint = request.get("routes_fingerprint")
    constraints: dict[str, Any] = {}
    target_structure = None
    artifact_index: dict[tuple[str, str], dict[str, Any]] = {}
    inventory = None
    inventory_index: dict[str, dict[str, Any]] = {}
    inventory_findings: list[dict[str, Any]] = []
    try:
        if request.get("schema_version") != SCHEMA_VERSION or request.get(
            "workflow"
        ) != WORKFLOW:
            raise InputFailure("schema_version/workflow 不匹配。")
        if request.get("input_profile") not in INPUT_PROFILES:
            if str(request.get("input_profile")).lower() in {"pickle", "pkl"}:
                top_findings.append(
                    finding("E-PICKLE-INPUT-001", "error", "input_profile")
                )
            else:
                raise InputFailure("input_profile 不受控。")
        if (
            options.get("comparison_mode") != "dimensions_only"
            or options.get("preserve_backend_order") is not True
        ):
            raise InputFailure(
                "首版必须显式使用 dimensions_only 和 preserve_backend_order=true。"
            )
        if not valid_sha256(source.get("content_sha256")):
            top_findings.append(
                finding("E-INPUT-HASH-001", "error", "source.content_sha256")
            )
        if SECRET_RE.search(canonical_json(request)):
            raise InputFailure("请求中检测到疑似凭证。")
        normalized_routes, route_input_findings = normalize_routes(request)
        top_findings.extend(route_input_findings)
        if not valid_sha256(routes_fingerprint) or routes_fingerprint != sha256_json(
            request.get("routes")
        ):
            top_findings.append(
                finding("E-INPUT-HASH-001", "error", "routes_fingerprint")
            )
        constraints = normalize_constraints(request.get("constraints"), toolkit)
        target = request.get("target")
        if isinstance(target, dict):
            _, target_structure = parse_molecule(
                target.get("standardized_structure")
                or target.get("reported_structure"),
                toolkit,
            )
            if (
                target.get("standardized_structure")
                or target.get("reported_structure")
            ) and target_structure is None:
                top_findings.append(
                    finding(
                        "E-MOLECULE-STRUCTURE-001",
                        "error",
                        "target",
                    )
                )
        artifact_index, artifact_findings = build_artifact_index(
            request.get("step_artifacts")
        )
        top_findings.extend(artifact_findings)
        inventory, inventory_index, inventory_findings = normalize_inventory(
            request.get("inventory_snapshot"), toolkit
        )
        if not any(item["severity"] == "error" for item in top_findings):
            total_nodes = 0
            for route in normalized_routes:
                processed = process_route(
                    route,
                    target_structure=target_structure,
                    artifact_index=artifact_index,
                    inventory=inventory,
                    inventory_index=inventory_index,
                    inventory_findings=inventory_findings,
                    constraints=constraints,
                    source_license=source.get("license"),
                    toolkit=toolkit,
                )
                total_nodes += processed["node_count"]
                routes.append(processed)
            if total_nodes > MAX_TOTAL_NODES:
                top_findings.append(
                    finding(
                        "E-RESOURCE-LIMIT-001",
                        "error",
                        "routes",
                        evidence={"total_nodes": total_nodes},
                    )
                )
                for route in routes:
                    route["findings"].append(top_findings[-1])
                    finalize_route(route)
            actual_step_keys = {
                (route["route_id"], step["step_id"])
                for route in routes
                for step in route["step_reviews"]
            }
            for route_id, step_id in sorted(set(artifact_index) - actual_step_keys):
                top_findings.append(
                    finding(
                        "E-STEP-HASH-MISMATCH-001",
                        "error",
                        "step_artifacts",
                        evidence={"route_id": route_id, "step_id": step_id},
                    )
                )
    except ResourceFailure as error:
        top_findings.append(
            finding(
                "E-RESOURCE-LIMIT-001",
                "error",
                "$",
                evidence=str(error),
            )
        )
    except (InputFailure, DependencyFailure) as error:
        top_findings.append(
            finding(
                "E-INPUT-SCHEMA-001",
                "error",
                "$",
                evidence=str(error),
            )
        )
    duplicates = duplicate_groups(routes)
    all_findings = list(top_findings)
    for route in routes:
        all_findings.extend(
            {"route_id": route["route_id"], **item}
            for item in route["findings"]
        )
    errors = [item for item in all_findings if item["severity"] == "error"]
    warnings = [item for item in all_findings if item["severity"] == "warning"]
    review_queue = []
    for route in routes:
        if route["disposition"] != "ready_for_expert_review":
            review_queue.append(
                {
                    "route_id": route["route_id"],
                    "step_id": None,
                    "reason_codes": route["human_review_required"]
                    or sorted(
                        {
                            item["code"]
                            for item in route["findings"]
                            if item["severity"] == "error"
                        }
                    ),
                }
            )
        for step in route["step_reviews"]:
            if step["review_required"]:
                review_queue.append(
                    {
                        "route_id": route["route_id"],
                        "step_id": step["step_id"],
                        "reason_codes": step["review_required"],
                    }
                )
    disposition_counts = {
        value: sum(route["disposition"] == value for route in routes)
        for value in sorted(DISPOSITIONS)
    }
    target_assessment = {
        "reported": (
            request.get("target") if isinstance(request.get("target"), dict) else None
        ),
        "canonical_structure": target_structure,
        "route_root_structures": sorted(
            {
                route["target_structure"]
                for route in routes
                if route["target_structure"]
            }
        ),
    }
    comparison_dimensions = [
        {
            "route_id": route["route_id"],
            "backend_rank": route["backend_metadata"]["backend_rank"],
            "backend_score": route["backend_metadata"]["backend_score"],
            "topology_status": route["topology_status"],
            "step_count": route["step_count"],
            "longest_linear_sequence": route["longest_linear_sequence"],
            "branch_count": route["branch_count"],
            "terminal_precursor_count": len(route["terminal_precursors"]),
            "inventory_coverage": route["inventory_coverage"],
            "precedent_coverage_by_level": route[
                "precedent_coverage_by_level"
            ],
            "exact_or_transformation_coverage": route[
                "exact_or_transformation_coverage"
            ],
            "weakest_step_count": len(route["weakest_steps"]),
            "disposition": route["disposition"],
        }
        for route in routes
    ]
    document = {
        "schema_version": SCHEMA_VERSION,
        "workflow": WORKFLOW,
        "ruleset_version": RULESET_VERSION,
        "generated_at_utc": generated_at_utc or now_utc(),
        "tool_versions": tool_versions(toolkit),
        "source_record": public_source,
        "routes_fingerprint": routes_fingerprint,
        "options": normalized_options,
        "constraints": constraints,
        "target_assessment": target_assessment,
        "input_summary": {
            "input_routes": len(request.get("routes") or [])
            if isinstance(request.get("routes"), list)
            else 0,
            "output_routes": len(routes),
            "total_nodes": sum(route["node_count"] for route in routes),
            "total_steps": sum(route["step_count"] for route in routes),
            "disposition_counts": disposition_counts,
            "record_count_conserved": (
                isinstance(request.get("routes"), list)
                and len(request["routes"]) == len(routes)
            ),
        },
        "route_summaries": routes,
        "duplicate_route_groups": duplicates,
        "comparison_dimensions": comparison_dimensions,
        "review_queue": review_queue,
        "errors": errors,
        "warnings": warnings,
        "notices": [
            "本输出是路线证据评审，不是可行性、安全、最优性或实验执行批准。",
            "backend score、库存声明、相似先例、条件和产率按各自来源保留，不生成默认综合总分。",
        ],
        "runtime_seconds": round(time.perf_counter() - started, 6),
    }
    document["result_fingerprint"] = stable_document_fingerprint(document)
    return document


def read_request(path: Path) -> dict[str, Any]:
    if path.suffix.lower() in {".pkl", ".pickle"}:
        raise InputFailure(RULE_MESSAGES["E-PICKLE-INPUT-001"])
    raw = path.read_text(encoding="utf-8")
    if SECRET_RE.search(raw):
        raise InputFailure("输入文件中检测到疑似凭证。")
    value = json.loads(raw)
    if not isinstance(value, dict):
        raise InputFailure("输入顶层必须是 JSON object。")
    return value


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args(argv)
    try:
        request = read_request(args.input)
        document = process_request(request)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(
            json.dumps(document, ensure_ascii=False, indent=2) + "\n",
            encoding="utf-8",
        )
        return 0 if not document["errors"] and all(
            route["disposition"] != "blocked"
            for route in document["route_summaries"]
        ) else 1
    except Exception as error:
        print(f"review-routes failed: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
