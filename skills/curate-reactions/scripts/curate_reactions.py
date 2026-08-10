#!/usr/bin/env python3
"""Non-destructive curation and quality review for structured reactions."""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.metadata
import json
import re
import sys
from collections import Counter, defaultdict
from collections.abc import Sequence
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

SCHEMA_VERSION = "1.0.0"
WORKFLOW = "curate-reactions"
RULESET_VERSION = "1.0.0"
MAX_RECORDS = 5000
MAX_RECORD_BYTES = 2 * 1024 * 1024
MAX_INPUT_BYTES = 100 * 1024 * 1024

CURATION_STATUSES = {"completed", "partial", "not_run", "error"}
DISPOSITIONS = {"ready_for_search", "review_required", "rejected"}
INPUT_PROFILES = {
    "ord_dataset",
    "ord_reaction",
    "reaction_smiles",
    "tabular",
}
PARTICIPANT_ROLES = {
    "reactant",
    "reagent",
    "catalyst",
    "solvent",
    "internal_standard",
    "product",
    "unknown",
}
TEMPORAL_TOP_LEVEL_KEYS = {
    "generated_at_utc",
    "runtime_seconds",
    "result_fingerprint",
}
SECRET_RE = re.compile(
    r"ark-[A-Za-z0-9_-]{12,}|"
    r"Bearer\s+[A-Za-z0-9._~+/=-]{12,}|"
    r"(?:Authorization|Cookie|Token|Api[_ -]?Key)\s*[:=]\s*\S{12,}",
    re.IGNORECASE,
)

RULE_MESSAGES = {
    "E-INPUT-SCHEMA-001": "输入顶层字段或枚举不符合冻结合同。",
    "E-INPUT-HASH-001": "来源缺少有效 SHA-256。",
    "E-RECORD-ID-001": "record_id 缺失或批内重复。",
    "E-ORD-PARSE-001": "ORD 记录无法按固定 Schema 解析。",
    "E-ORD-VALIDATION-001": "ORD 官方校验返回阻断错误。",
    "E-REACTION-SIDES-001": "反应缺少可识别的输入或输出侧。",
    "E-REACTION-SMILES-001": "reaction SMILES 不是可解析的两段或三段形式。",
    "E-UPSTREAM-FINGERPRINT-001": "上游 artifact 指纹缺失或与内容不匹配。",
    "E-RESOURCE-LIMIT-001": "输入超过首版资源上限。",
    "W-PARTICIPANT-STRUCTURE-001": "参与物结构不可解析或上游状态拒绝。",
    "W-PARTICIPANT-FORM-001": "报告形式、标准化形式或 parent 存在差异。",
    "W-ROLE-CONFLICT-001": "来源角色与可确定的参与性诊断冲突。",
    "W-ROLE-UNKNOWN-001": "参与物角色不能确定。",
    "W-YIELD-RANGE-001": "百分比产率超出 0 到 100。",
    "W-YIELD-FRACTION-001": "百分比产率疑似混用 0 到 1 小数。",
    "W-YIELD-CONFLICT-001": "同一产物存在冲突的产率记录。",
    "W-ANALYSIS-LINK-001": "产物测量没有关联分析记录。",
    "W-DUPLICATE-EXACT-001": "记录与其他记录在 exact_record 视图重复。",
    "W-DUPLICATE-TRANSFORMATION-001": "报告形式的反应转化重复。",
    "H-DUPLICATE-PARENT-001": "记录只在 parent 转化视图相同，禁止自动合并。",
    "W-BALANCE-ATOM-001": "显式反应结构的逐元素计数不平衡。",
    "W-BALANCE-CHARGE-001": "显式反应结构的形式电荷不平衡。",
    "H-BALANCE-INCOMPLETE-001": "守恒检查缺少计量数、共反应物或副产物前提。",
    "W-MAPPING-FAILED-001": "请求了原子映射，但首版核心未运行或 adapter 失败。",
    "H-MAPPING-LOW-CONFIDENCE-001": "原子映射置信度不足或位于适用域外。",
    "W-PROCESS-MISSING-001": "用途要求的条件、装置、观察或后处理字段缺失。",
    "W-ORD-UNCLASSIFIED-001": "ORD 返回未分类 warning，原文已保留。",
    "W-REACTION-NO-CHANGE-001": "输入和输出结构集合相同，需确认是否为有效反应记录。",
    "H-UPSTREAM-REVIEW-001": "参与物继承上游人工复核状态。",
}


class InputFailure(RuntimeError):
    """Raised for a top-level input failure."""


def now_utc() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def canonical_json(value: Any) -> str:
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"))


def sha256_json(value: Any) -> str:
    return hashlib.sha256(canonical_json(value).encode("utf-8")).hexdigest()


def stable_document_fingerprint(document: dict[str, Any]) -> str:
    payload = {
        key: value
        for key, value in document.items()
        if key not in TEMPORAL_TOP_LEVEL_KEYS
    }
    return sha256_json(payload)


def upstream_fingerprint(document: dict[str, Any]) -> str:
    payload = {
        key: value
        for key, value in document.items()
        if key not in {"generated_at_utc", "result_fingerprint"}
    }
    return sha256_json(payload)


def finding(
    code: str,
    severity: str,
    field_path: str,
    *,
    detail: str | None = None,
    raw_message: str | None = None,
    evidence: Sequence[dict[str, Any]] | None = None,
) -> dict[str, Any]:
    item: dict[str, Any] = {
        "code": code,
        "severity": severity,
        "field_path": field_path,
        "message": RULE_MESSAGES[code],
        "evidence": list(evidence or []),
    }
    if detail is not None:
        item["detail"] = detail
    if raw_message is not None:
        item["raw_message"] = raw_message
    return item


def load_toolkit() -> dict[str, Any]:
    try:
        import rdkit
        from rdkit import Chem
    except ImportError as exc:
        raise InputFailure(
            "缺少 rdkit==2025.9.2；请在隔离环境安装 scripts/requirements.txt"
        ) from exc
    rdkit_version = importlib.metadata.version("rdkit")
    if rdkit_version != "2025.9.2":
        raise InputFailure(f"rdkit 版本必须为 2025.9.2，当前为 {rdkit_version}")
    try:
        from google.protobuf.json_format import MessageToDict, ParseDict
        from ord_schema import message_helpers, validations
        from ord_schema.proto import dataset_pb2, reaction_pb2
    except ImportError as exc:
        raise InputFailure(
            "缺少 ord-schema==0.8.3；请在隔离环境安装 scripts/requirements.txt"
        ) from exc
    ord_version = importlib.metadata.version("ord-schema")
    if ord_version != "0.8.3":
        raise InputFailure(f"ord-schema 版本必须为 0.8.3，当前为 {ord_version}")
    return {
        "rdkit": rdkit,
        "Chem": Chem,
        "MessageToDict": MessageToDict,
        "ParseDict": ParseDict,
        "message_helpers": message_helpers,
        "validations": validations,
        "dataset_pb2": dataset_pb2,
        "reaction_pb2": reaction_pb2,
        "rdkit_version": rdkit_version,
        "ord_version": ord_version,
    }


def validate_sha256(value: Any) -> bool:
    return isinstance(value, str) and bool(re.fullmatch(r"[0-9a-fA-F]{64}", value))


def split_reaction_smiles(value: Any) -> tuple[list[str], list[str], list[str]]:
    if not isinstance(value, str) or not value.strip():
        raise ValueError("reaction SMILES 为空")
    text = value.strip()
    if text.count(">") != 2:
        raise ValueError("reaction SMILES 必须包含两个 >")
    left, middle, right = text.split(">")
    inputs = [item for item in left.split(".") if item]
    agents = [item for item in middle.split(".") if item]
    outputs = [item for item in right.split(".") if item]
    if not inputs or not outputs:
        raise ValueError("反应必须同时包含输入和输出")
    return inputs, agents, outputs


def canonicalize_smiles(value: str, toolkit: dict[str, Any]) -> tuple[str | None, Any]:
    Chem = toolkit["Chem"]
    with toolkit["rdkit"].rdBase.BlockLogs():
        mol = Chem.MolFromSmiles(value)
    if mol is None:
        return None, None
    return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True), mol


def molecule_composition(mol: Any) -> tuple[Counter[str], int]:
    counts: Counter[str] = Counter()
    charge = 0
    for atom in mol.GetAtoms():
        counts[atom.GetSymbol()] += 1
        counts["H"] += int(atom.GetTotalNumHs(includeNeighbors=True))
        charge += int(atom.GetFormalCharge())
    if not counts["H"]:
        counts.pop("H", None)
    return counts, charge


def composition_delta(
    input_mols: Sequence[Any], output_mols: Sequence[Any]
) -> tuple[dict[str, int], int]:
    inputs: Counter[str] = Counter()
    outputs: Counter[str] = Counter()
    input_charge = 0
    output_charge = 0
    for mol in input_mols:
        counts, charge = molecule_composition(mol)
        inputs.update(counts)
        input_charge += charge
    for mol in output_mols:
        counts, charge = molecule_composition(mol)
        outputs.update(counts)
        output_charge += charge
    elements = sorted(set(inputs) | set(outputs))
    delta = {
        element: outputs[element] - inputs[element]
        for element in elements
        if outputs[element] != inputs[element]
    }
    return delta, output_charge - input_charge


def classify_ord_warning(message: str) -> str:
    if "analysis_key" in message and "Product measurements" in message:
        return "W-ANALYSIS-LINK-001"
    if "outside the expected range (0-100)" in message:
        return "W-YIELD-RANGE-001"
    if "Percentage values are 0-100, not fractions" in message:
        return "W-YIELD-FRACTION-001"
    return "W-ORD-UNCLASSIFIED-001"


def parse_ord_record(
    raw: dict[str, Any], toolkit: dict[str, Any]
) -> tuple[dict[str, Any] | None, str | None, list[dict[str, Any]]]:
    findings: list[dict[str, Any]] = []
    reaction = toolkit["reaction_pb2"].Reaction()
    try:
        toolkit["ParseDict"](raw, reaction, ignore_unknown_fields=False)
    except Exception as exc:  # noqa: BLE001 - protobuf exposes mixed parse errors
        findings.append(
            finding(
                "E-ORD-PARSE-001",
                "error",
                "ord_record",
                detail=f"{type(exc).__name__}: {exc}",
            )
        )
        return None, None, findings
    output = toolkit["validations"].validate_message(
        reaction, raise_on_error=False
    )
    for message in output.errors:
        findings.append(
            finding(
                "E-ORD-VALIDATION-001",
                "error",
                "ord_record",
                raw_message=str(message),
            )
        )
    for message in output.warnings:
        findings.append(
            finding(
                classify_ord_warning(str(message)),
                "warning",
                "ord_record",
                raw_message=str(message),
            )
        )
    try:
        reaction_smiles = toolkit["message_helpers"].get_reaction_smiles(
            reaction, generate_if_missing=True
        )
    except Exception as exc:  # noqa: BLE001 - ORD helper may wrap RDKit errors
        reaction_smiles = None
        findings.append(
            finding(
                "E-REACTION-SMILES-001",
                "error",
                "ord_record.identifiers",
                detail=f"无法从 ORD 生成 reaction SMILES：{exc}",
            )
        )
    normalized = toolkit["MessageToDict"](
        reaction,
        preserving_proto_field_name=True,
        use_integers_for_enums=False,
    )
    return normalized, reaction_smiles, findings


def extract_ord_yields(ord_record: Any) -> list[dict[str, Any]]:
    if not isinstance(ord_record, dict):
        return []
    yields = []
    for outcome_index, outcome in enumerate(ord_record.get("outcomes") or []):
        if not isinstance(outcome, dict):
            continue
        for product_index, product in enumerate(outcome.get("products") or []):
            if not isinstance(product, dict):
                continue
            product_id = f"outcome-{outcome_index + 1}-product-{product_index + 1}"
            identifiers = product.get("identifiers")
            if isinstance(identifiers, list):
                for identifier in identifiers:
                    if (
                        isinstance(identifier, dict)
                        and identifier.get("type") == "SMILES"
                        and identifier.get("value")
                    ):
                        product_id = str(identifier["value"])
                        break
            for measurement in product.get("measurements") or []:
                if (
                    not isinstance(measurement, dict)
                    or measurement.get("type") != "YIELD"
                ):
                    continue
                percentage = measurement.get("percentage")
                value = (
                    percentage.get("value")
                    if isinstance(percentage, dict)
                    else None
                )
                yields.append(
                    {
                        "value": value,
                        "units": "PERCENT",
                        "type": "reported",
                        "product_id": product_id,
                        "analysis_key": measurement.get("analysis_key"),
                        # ORD's official validator already reports missing links.
                        "analysis_required": False,
                    }
                )
    return yields


def build_upstream_index(
    artifacts: Sequence[Any],
) -> tuple[dict[str, dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    index: dict[str, dict[str, Any]] = {}
    metadata: list[dict[str, Any]] = []
    errors: list[dict[str, Any]] = []
    for position, value in enumerate(artifacts):
        artifact = value.get("artifact") if isinstance(value, dict) else None
        if artifact is None and isinstance(value, dict):
            artifact = value
        if not isinstance(artifact, dict):
            errors.append(
                finding(
                    "E-UPSTREAM-FINGERPRINT-001",
                    "error",
                    f"upstream_artifacts[{position}]",
                    detail="artifact 必须是 object",
                )
            )
            continue
        actual = artifact.get("result_fingerprint")
        expected = upstream_fingerprint(artifact)
        if not validate_sha256(actual) or actual != expected:
            errors.append(
                finding(
                    "E-UPSTREAM-FINGERPRINT-001",
                    "error",
                    f"upstream_artifacts[{position}].result_fingerprint",
                )
            )
            continue
        records = artifact.get("records")
        if not isinstance(records, list):
            errors.append(
                finding(
                    "E-UPSTREAM-FINGERPRINT-001",
                    "error",
                    f"upstream_artifacts[{position}].records",
                    detail="上游 artifact 缺少 records",
                )
            )
            continue
        for record in records:
            if isinstance(record, dict) and isinstance(record.get("id"), str):
                index[record["id"]] = record
        metadata.append(
            {
                "workflow": artifact.get("workflow"),
                "schema_version": artifact.get("schema_version"),
                "result_fingerprint": actual,
                "record_count": len(records),
            }
        )
    return index, metadata, errors


def source_participants(
    reaction_smiles: str,
) -> list[dict[str, Any]]:
    inputs, agents, outputs = split_reaction_smiles(reaction_smiles)
    participants = []
    for side, role, structures in (
        ("input", "reactant", inputs),
        ("input", "reagent", agents),
        ("output", "product", outputs),
    ):
        for index, structure in enumerate(structures):
            participants.append(
                {
                    "participant_id": f"{side}-{role}-{index + 1}",
                    "side": side,
                    "reported_role": role,
                    "original_structure": structure,
                }
            )
    return participants


def assess_participant(
    raw: dict[str, Any],
    upstream: dict[str, dict[str, Any]],
    toolkit: dict[str, Any],
    index: int,
) -> tuple[dict[str, Any], Any]:
    participant_id = raw.get("participant_id") or f"participant-{index + 1}"
    role = raw.get("reported_role", "unknown")
    side = raw.get("side")
    structure = raw.get("original_structure")
    upstream_id = raw.get("upstream_record_id")
    upstream_record = upstream.get(upstream_id) if isinstance(upstream_id, str) else None
    findings: list[dict[str, Any]] = []
    if role not in PARTICIPANT_ROLES:
        role = "unknown"
    if side not in {"input", "output"}:
        side = "output" if role == "product" else "input"
    standardized = None
    parent = None
    upstream_disposition = None
    upstream_review: list[Any] = []
    if upstream_record is not None:
        structure = structure or upstream_record.get("original_structure")
        standardized = upstream_record.get("standardized_structure")
        parent = upstream_record.get("parent_structure")
        upstream_disposition = upstream_record.get("disposition")
        upstream_review = list(upstream_record.get("human_review_required") or [])
        if upstream_disposition == "rejected":
            findings.append(
                finding(
                    "W-PARTICIPANT-STRUCTURE-001",
                    "warning",
                    f"participants[{index}].upstream_disposition",
                    detail="上游参与物已 rejected",
                )
            )
        if upstream_disposition == "review_required" or upstream_review:
            findings.append(
                finding(
                    "H-UPSTREAM-REVIEW-001",
                    "human_review",
                    f"participants[{index}]",
                    detail="; ".join(map(str, upstream_review)) or "上游要求人工复核",
                )
            )
    canonical, mol = (None, None)
    if isinstance(standardized, str):
        canonical, mol = canonicalize_smiles(standardized, toolkit)
    if canonical is None and isinstance(structure, str):
        canonical, mol = canonicalize_smiles(structure, toolkit)
    if canonical is None:
        findings.append(
            finding(
                "W-PARTICIPANT-STRUCTURE-001",
                "warning",
                f"participants[{index}].original_structure",
            )
        )
    if (
        isinstance(structure, str)
        and isinstance(standardized, str)
        and structure != standardized
    ) or (
        isinstance(parent, str)
        and isinstance(standardized, str)
        and parent != standardized
    ):
        findings.append(
            finding(
                "W-PARTICIPANT-FORM-001",
                "warning",
                f"participants[{index}]",
                detail="保留 reported/standardized/parent 差异",
            )
        )
    if role == "unknown":
        findings.append(
            finding(
                "W-ROLE-UNKNOWN-001",
                "warning",
                f"participants[{index}].reported_role",
            )
        )
    participation_status = (
        "product"
        if side == "output"
        else "not_assessed"
        if role in {"reagent", "catalyst", "solvent", "unknown"}
        else "contributes_product_atoms"
    )
    if side == "output" and role != "product":
        findings.append(
            finding(
                "W-ROLE-CONFLICT-001",
                "warning",
                f"participants[{index}].reported_role",
                detail="输出侧参与物未报告为 product",
            )
        )
    if side == "input" and role == "product":
        findings.append(
            finding(
                "W-ROLE-CONFLICT-001",
                "warning",
                f"participants[{index}].reported_role",
                detail="输入侧参与物报告为 product",
            )
        )
    assessment = {
        "participant_id": participant_id,
        "side": side,
        "reported_role": role,
        "reported_form": structure,
        "standardized_form": canonical,
        "parent_form": parent,
        "upstream_record_id": upstream_id,
        "upstream_disposition": upstream_disposition,
        "upstream_human_review_required": upstream_review,
        "participation_status": participation_status,
        "role_status": (
            "conflict"
            if any(item["code"] == "W-ROLE-CONFLICT-001" for item in findings)
            else "not_assessed"
            if role == "unknown"
            else "consistent"
        ),
        "findings": findings,
    }
    return assessment, mol


def assess_yields(raw: dict[str, Any]) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    yields = raw.get("yields")
    if yields is None and raw.get("yield_percent") is not None:
        yields = [
            {
                "value": raw.get("yield_percent"),
                "units": "PERCENT",
                "type": "reported",
            }
        ]
    if not isinstance(yields, list):
        yields = []
    normalized = []
    findings: list[dict[str, Any]] = []
    by_product: defaultdict[str, set[float]] = defaultdict(set)
    for index, entry in enumerate(yields):
        if not isinstance(entry, dict):
            continue
        value = entry.get("value")
        item = {
            "value": value,
            "units": entry.get("units", "PERCENT"),
            "type": entry.get("type", "reported"),
            "product_id": entry.get("product_id"),
            "analysis_key": entry.get("analysis_key"),
        }
        normalized.append(item)
        if isinstance(value, (int, float)) and not isinstance(value, bool):
            if value < 0 or value > 100:
                findings.append(
                    finding(
                        "W-YIELD-RANGE-001",
                        "warning",
                        f"yields[{index}].value",
                        detail=str(value),
                    )
                )
            elif 0 < value < 1:
                findings.append(
                    finding(
                        "W-YIELD-FRACTION-001",
                        "warning",
                        f"yields[{index}].value",
                        detail=str(value),
                    )
                )
            by_product[str(entry.get("product_id"))].add(float(value))
        if entry.get("analysis_required") is True and not entry.get("analysis_key"):
            findings.append(
                finding(
                    "W-ANALYSIS-LINK-001",
                    "warning",
                    f"yields[{index}].analysis_key",
                )
            )
    for product_id, values in by_product.items():
        if len(values) > 1:
            findings.append(
                finding(
                    "W-YIELD-CONFLICT-001",
                    "warning",
                    "yields",
                    detail=f"product_id={product_id}, values={sorted(values)}",
                )
            )
    return findings, {"measurements": normalized, "status": "completed"}


def assess_record(
    raw: dict[str, Any],
    upstream: dict[str, dict[str, Any]],
    toolkit: dict[str, Any],
) -> dict[str, Any]:
    record_id = raw.get("record_id")
    original_hash = sha256_json(raw)
    findings: list[dict[str, Any]] = []
    ord_record = raw.get("ord_record")
    reaction_smiles = raw.get("reaction_smiles")
    if isinstance(reaction_smiles, dict):
        reaction_smiles = reaction_smiles.get("reported")
    normalized_ord = None
    if isinstance(ord_record, dict):
        normalized_ord, generated_smiles, ord_findings = parse_ord_record(
            ord_record, toolkit
        )
        findings.extend(ord_findings)
        reaction_smiles = reaction_smiles or generated_smiles
    try:
        input_smiles, agent_smiles, output_smiles = split_reaction_smiles(
            reaction_smiles
        )
    except ValueError as exc:
        findings.append(
            finding(
                "E-REACTION-SMILES-001",
                "error",
                "reaction_smiles",
                detail=str(exc),
            )
        )
        input_smiles, agent_smiles, output_smiles = [], [], []
    if not input_smiles or not output_smiles:
        findings.append(
            finding(
                "E-REACTION-SIDES-001",
                "error",
                "reaction_smiles",
            )
        )
    participants = raw.get("participants")
    if not isinstance(participants, list) or not participants:
        try:
            participants = source_participants(str(reaction_smiles))
        except ValueError:
            participants = []
    assessments = []
    input_mols = []
    output_mols = []
    for index, participant in enumerate(participants):
        if not isinstance(participant, dict):
            participant = {}
        assessment, mol = assess_participant(
            participant, upstream, toolkit, index
        )
        assessments.append(assessment)
        findings.extend(assessment["findings"])
        if mol is not None:
            if assessment["side"] == "output":
                output_mols.append(mol)
            elif assessment["reported_role"] == "reactant":
                input_mols.append(mol)
    yield_source = raw
    if (
        raw.get("yields") is None
        and raw.get("yield_percent") is None
        and normalized_ord is not None
    ):
        yield_source = dict(raw)
        yield_source["yields"] = extract_ord_yields(normalized_ord)
    yield_findings, yield_assessment = assess_yields(yield_source)
    findings.extend(yield_findings)
    balance_assessment: dict[str, Any] = {
        "status": "not_assessed",
        "assumption": "each_listed_reactant_once",
        "element_delta": {},
        "formal_charge_delta": 0,
    }
    if input_mols and output_mols:
        element_delta, charge_delta = composition_delta(input_mols, output_mols)
        balance_assessment.update(
            {
                "status": "completed",
                "element_delta": element_delta,
                "formal_charge_delta": charge_delta,
            }
        )
        if element_delta:
            findings.append(
                finding(
                    "W-BALANCE-ATOM-001",
                    "warning",
                    "balance_assessment.element_delta",
                    detail=canonical_json(element_delta),
                )
            )
        if charge_delta:
            findings.append(
                finding(
                    "W-BALANCE-CHARGE-001",
                    "warning",
                    "balance_assessment.formal_charge_delta",
                    detail=str(charge_delta),
                )
            )
        if raw.get("stoichiometry_complete") is not True:
            findings.append(
                finding(
                    "H-BALANCE-INCOMPLETE-001",
                    "human_review",
                    "balance_assessment",
                )
            )
    mapping_assessment = {
        "requested": False,
        "status": "not_run",
        "backend": None,
        "confidence": None,
    }
    process = raw.get("process")
    if isinstance(process, dict) and process.get("required") is True:
        missing = [
            key
            for key in ("conditions", "setup", "observations", "workups")
            if not process.get(key)
        ]
        if missing:
            findings.append(
                finding(
                    "W-PROCESS-MISSING-001",
                    "warning",
                    "process",
                    detail=", ".join(missing),
                )
            )
    canonical_inputs = []
    canonical_agents = []
    canonical_outputs = []
    for values, target in (
        (input_smiles, canonical_inputs),
        (agent_smiles, canonical_agents),
        (output_smiles, canonical_outputs),
    ):
        for value in values:
            canonical, _ = canonicalize_smiles(value, toolkit)
            if canonical:
                target.append(canonical)
    canonical_reaction = (
        ".".join(sorted(canonical_inputs))
        + ">"
        + ".".join(sorted(canonical_agents))
        + ">"
        + ".".join(sorted(canonical_outputs))
        if canonical_inputs and canonical_outputs
        else None
    )
    if (
        canonical_inputs
        and canonical_outputs
        and sorted(canonical_inputs) == sorted(canonical_outputs)
    ):
        findings.append(
            finding(
                "W-REACTION-NO-CHANGE-001",
                "warning",
                "reaction_smiles",
            )
        )
    severities = {item["severity"] for item in findings}
    if "error" in severities:
        curation_status = "error"
        disposition = "rejected"
    elif findings:
        curation_status = "partial"
        disposition = "review_required"
    else:
        curation_status = "completed"
        disposition = "ready_for_search"
    exact_payload = {
        "reaction_smiles": canonical_reaction,
        "participants": [
            {
                "side": item["side"],
                "reported_role": item["reported_role"],
                "standardized_form": item["standardized_form"],
            }
            for item in assessments
        ],
        "yields": yield_assessment["measurements"],
        "ord_record": normalized_ord,
    }
    parent_inputs = sorted(
        item["parent_form"]
        for item in assessments
        if item["side"] == "input" and isinstance(item["parent_form"], str)
    )
    parent_outputs = sorted(
        item["parent_form"]
        for item in assessments
        if item["side"] == "output" and isinstance(item["parent_form"], str)
    )
    parent_key = (
        ".".join(parent_inputs) + ">>" + ".".join(parent_outputs)
        if parent_inputs and parent_outputs
        else None
    )
    return {
        "record_id": record_id,
        "source_locator": raw.get("source_locator"),
        "original_record_hash": original_hash,
        "ord_record": normalized_ord,
        "reaction_smiles": {
            "reported": reaction_smiles,
            "canonical_unmapped": canonical_reaction,
        },
        "participant_assessments": assessments,
        "role_assessment": {
            "status": (
                "conflict"
                if any(item["role_status"] == "conflict" for item in assessments)
                else "review_required"
                if any(
                    item["role_status"] in {"ambiguous", "not_assessed"}
                    for item in assessments
                )
                else "consistent"
            )
        },
        "yield_assessment": yield_assessment,
        "balance_assessment": balance_assessment,
        "mapping_assessment": mapping_assessment,
        "duplicate_memberships": [],
        "curation_status": curation_status,
        "findings": findings,
        "disposition": disposition,
        "human_review_required": sorted(
            {
                item["code"]
                for item in findings
                if item["severity"] == "human_review"
            }
        ),
        "_duplicate_keys": {
            "exact_record": sha256_json(exact_payload),
            "reported_transformation": canonical_reaction,
            "parent_transformation_candidate": parent_key,
        },
    }


def apply_duplicate_groups(
    records: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    groups_by_view: dict[str, defaultdict[str, list[dict[str, Any]]]] = {
        view: defaultdict(list)
        for view in (
            "exact_record",
            "reported_transformation",
            "parent_transformation_candidate",
        )
    }
    for record in records:
        for view, key in record["_duplicate_keys"].items():
            if key:
                groups_by_view[view][key].append(record)
    groups = []
    code_for_view = {
        "exact_record": "W-DUPLICATE-EXACT-001",
        "reported_transformation": "W-DUPLICATE-TRANSFORMATION-001",
        "parent_transformation_candidate": "H-DUPLICATE-PARENT-001",
    }
    for view in (
        "exact_record",
        "reported_transformation",
        "parent_transformation_candidate",
    ):
        for key in sorted(groups_by_view[view]):
            members = groups_by_view[view][key]
            if len(members) < 2:
                continue
            group_id = f"{view}:{hashlib.sha256(key.encode()).hexdigest()[:16]}"
            member_ids = [item["record_id"] for item in members]
            groups.append(
                {
                    "group_id": group_id,
                    "view": view,
                    "record_ids": member_ids,
                    "automatic_action": "none",
                }
            )
            code = code_for_view[view]
            severity = (
                "human_review"
                if view == "parent_transformation_candidate"
                else "warning"
            )
            for record in members:
                record["duplicate_memberships"].append(group_id)
                record["findings"].append(
                    finding(
                        code,
                        severity,
                        "duplicate_memberships",
                        detail=f"{view}: {member_ids}",
                    )
                )
                if record["disposition"] == "ready_for_search":
                    record["disposition"] = "review_required"
                    record["curation_status"] = "partial"
                if severity == "human_review":
                    record["human_review_required"] = sorted(
                        set(record["human_review_required"]) | {code}
                    )
    for record in records:
        record.pop("_duplicate_keys", None)
    return groups


def process_request(
    request: dict[str, Any],
    *,
    generated_at_utc: str | None = None,
) -> dict[str, Any]:
    toolkit = load_toolkit()
    top_errors: list[dict[str, Any]] = []
    if request.get("schema_version") != SCHEMA_VERSION:
        top_errors.append(
            finding("E-INPUT-SCHEMA-001", "error", "schema_version")
        )
    if request.get("workflow") != WORKFLOW:
        top_errors.append(finding("E-INPUT-SCHEMA-001", "error", "workflow"))
    if request.get("input_profile") not in INPUT_PROFILES:
        top_errors.append(
            finding("E-INPUT-SCHEMA-001", "error", "input_profile")
        )
    source = request.get("source")
    if not isinstance(source, dict) or not validate_sha256(
        source.get("content_sha256")
    ):
        top_errors.append(
            finding("E-INPUT-HASH-001", "error", "source.content_sha256")
        )
    options = request.get("options")
    if not isinstance(options, dict):
        options = {}
    normalized_options = {
        "participant_view": options.get("participant_view", "reported_form"),
        "atom_mapping": options.get("atom_mapping", "off"),
        "balance_check": options.get("balance_check", "diagnostic"),
        "duplicate_views": options.get(
            "duplicate_views",
            [
                "exact_record",
                "reported_transformation",
                "parent_transformation_candidate",
            ],
        ),
        "preserve_original": True,
        "network_access": False,
        "automatic_writeback": False,
    }
    if (
        normalized_options["participant_view"] != "reported_form"
        or normalized_options["atom_mapping"] != "off"
        or normalized_options["balance_check"] != "diagnostic"
    ):
        top_errors.append(
            finding(
                "E-INPUT-SCHEMA-001",
                "error",
                "options",
                detail="首版只允许 reported_form / atom_mapping=off / diagnostic",
            )
        )
        normalized_options["participant_view"] = "reported_form"
        normalized_options["atom_mapping"] = "off"
        normalized_options["balance_check"] = "diagnostic"
    records = request.get("records")
    if not isinstance(records, list):
        records = []
        top_errors.append(finding("E-INPUT-SCHEMA-001", "error", "records"))
    if len(records) > MAX_RECORDS:
        top_errors.append(
            finding("E-RESOURCE-LIMIT-001", "error", "records")
        )
    run_fatal = any(
        item["code"]
        in {
            "E-INPUT-SCHEMA-001",
            "E-INPUT-HASH-001",
            "E-RESOURCE-LIMIT-001",
        }
        for item in top_errors
    )
    seen_ids: set[str] = set()
    invalid_record_ids: set[int] = set()
    for index, record in enumerate(records):
        if len(canonical_json(record).encode("utf-8")) > MAX_RECORD_BYTES:
            top_errors.append(
                finding(
                    "E-RESOURCE-LIMIT-001",
                    "error",
                    f"records[{index}]",
                )
            )
            invalid_record_ids.add(index)
        record_id = record.get("record_id") if isinstance(record, dict) else None
        if (
            not isinstance(record_id, str)
            or not record_id
            or record_id in seen_ids
        ):
            top_errors.append(
                finding(
                    "E-RECORD-ID-001",
                    "error",
                    f"records[{index}].record_id",
                )
            )
            invalid_record_ids.add(index)
        else:
            seen_ids.add(record_id)
    upstream, upstream_metadata, upstream_errors = build_upstream_index(
        request.get("upstream_artifacts") or []
    )
    top_errors.extend(upstream_errors)
    processed = []
    records_to_process = [] if len(records) > MAX_RECORDS else records
    if run_fatal:
        invalid_record_ids.update(range(len(records_to_process)))
    for index, raw in enumerate(records_to_process):
        if not isinstance(raw, dict):
            continue
        if index in invalid_record_ids:
            processed.append(
                {
                    "record_id": f"__invalid_record_{index + 1}",
                    "source_locator": {
                        "original_record_id": raw.get("record_id"),
                        "source_locator": raw.get("source_locator"),
                    },
                    "original_record_hash": sha256_json(raw),
                    "ord_record": None,
                    "reaction_smiles": {
                        "reported": raw.get("reaction_smiles"),
                        "canonical_unmapped": None,
                    },
                    "participant_assessments": [],
                    "role_assessment": {"status": "not_assessed"},
                    "yield_assessment": {
                        "measurements": [],
                        "status": "not_run",
                    },
                    "balance_assessment": {
                        "status": "not_assessed",
                        "assumption": "none",
                        "element_delta": {},
                        "formal_charge_delta": 0,
                    },
                    "mapping_assessment": {
                        "requested": False,
                        "status": "not_run",
                        "backend": None,
                        "confidence": None,
                    },
                    "duplicate_memberships": [],
                    "curation_status": "error",
                    "findings": [
                        item
                        for item in top_errors
                        if run_fatal
                        or item["field_path"].startswith(f"records[{index}]")
                    ],
                    "disposition": "rejected",
                    "human_review_required": [],
                    "_duplicate_keys": {
                        "exact_record": None,
                        "reported_transformation": None,
                        "parent_transformation_candidate": None,
                    },
                }
            )
            continue
        processed.append(assess_record(raw, upstream, toolkit))
    duplicate_groups = apply_duplicate_groups(processed)
    all_findings = list(top_errors)
    for record in processed:
        all_findings.extend(
            {"record_id": record["record_id"], **item}
            for item in record["findings"]
        )
    errors = [item for item in all_findings if item["severity"] == "error"]
    warnings = [item for item in all_findings if item["severity"] == "warning"]
    human = [
        item for item in all_findings if item["severity"] == "human_review"
    ]
    review_queue = [
        {
            "record_id": record["record_id"],
            "required_action": "human_review",
            "reason_codes": sorted(
                {
                    item["code"]
                    for item in record["findings"]
                    if item["severity"] in {"warning", "human_review"}
                }
            ),
        }
        for record in processed
        if record["disposition"] == "review_required"
    ]
    dispositions = {
        status: sum(record["disposition"] == status for record in processed)
        for status in sorted(DISPOSITIONS)
    }
    statuses = {
        status: sum(record["curation_status"] == status for record in processed)
        for status in sorted(CURATION_STATUSES)
    }
    document = {
        "schema_version": SCHEMA_VERSION,
        "workflow": WORKFLOW,
        "ruleset_version": RULESET_VERSION,
        "generated_at_utc": generated_at_utc or now_utc(),
        "tool_versions": {
            "rdkit": toolkit["rdkit_version"],
            "ord-schema": toolkit["ord_version"],
        },
        "options": normalized_options,
        "source_record": {
            "identifier": source.get("identifier") if isinstance(source, dict) else None,
            "content_sha256": (
                source.get("content_sha256") if isinstance(source, dict) else None
            ),
            "license": source.get("license") if isinstance(source, dict) else None,
        },
        "upstream_artifacts": upstream_metadata,
        "input_summary": {
            "total_records": len(records),
            "output_records": len(processed),
            "disposition_counts": dispositions,
            "curation_status_counts": statuses,
        },
        "records": processed,
        "duplicate_groups": duplicate_groups,
        "review_queue": review_queue,
        "errors": errors,
        "warnings": warnings,
        "notices": [
            "ready_for_search 只表示通过当前结构化检索数据门槛，不表示适合建模。",
            "元素/电荷平衡和原子映射诊断不证明反应正确、可复现、安全或可执行。",
            "本工作流不删除、合并、覆盖或写回原始反应记录。",
        ],
        "human_review_required": human,
    }
    document["result_fingerprint"] = stable_document_fingerprint(document)
    serialized = canonical_json(document)
    if SECRET_RE.search(serialized):
        raise InputFailure("输出中检测到疑似凭证，已停止写出")
    return document


def load_request(path: Path, toolkit: dict[str, Any]) -> dict[str, Any]:
    if path.stat().st_size > MAX_INPUT_BYTES:
        raise InputFailure("输入文件超过 100 MiB")
    content_hash = hashlib.sha256(path.read_bytes()).hexdigest()
    suffixes = path.suffixes
    if path.suffix.lower() == ".csv":
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            rows = list(csv.DictReader(handle))
        if not rows or not {"record_id", "reaction_smiles"}.issubset(
            rows[0].keys()
        ):
            raise InputFailure("CSV 必须包含 record_id,reaction_smiles")
        records = []
        for row in rows:
            record: dict[str, Any] = {
                "record_id": row.get("record_id"),
                "reaction_smiles": row.get("reaction_smiles"),
            }
            if row.get("yield_percent"):
                try:
                    record["yield_percent"] = float(row["yield_percent"])
                except ValueError:
                    record["yield_percent"] = row["yield_percent"]
            records.append(record)
        return {
            "schema_version": SCHEMA_VERSION,
            "workflow": WORKFLOW,
            "input_profile": "tabular",
            "source": {
                "identifier": path.name,
                "content_sha256": content_hash,
            },
            "options": {},
            "upstream_artifacts": [],
            "records": records,
        }
    if suffixes[-2:] == [".pb", ".gz"] or path.suffix.lower() == ".pb":
        dataset = toolkit["message_helpers"].load_message(
            str(path), toolkit["dataset_pb2"].Dataset
        )
        records = [
            {
                "record_id": reaction.reaction_id or f"reaction-{index + 1}",
                "ord_record": toolkit["MessageToDict"](
                    reaction,
                    preserving_proto_field_name=True,
                    use_integers_for_enums=False,
                ),
            }
            for index, reaction in enumerate(dataset.reactions)
        ]
        return {
            "schema_version": SCHEMA_VERSION,
            "workflow": WORKFLOW,
            "input_profile": "ord_dataset",
            "source": {
                "identifier": dataset.dataset_id or path.name,
                "content_sha256": content_hash,
            },
            "options": {},
            "upstream_artifacts": [],
            "records": records,
        }
    try:
        document = json.loads(path.read_text(encoding="utf-8"))
    except (UnicodeError, json.JSONDecodeError) as exc:
        raise InputFailure(f"无法读取 JSON：{exc}") from exc
    if not isinstance(document, dict):
        raise InputFailure("JSON 顶层必须是 object")
    if document.get("workflow") == WORKFLOW:
        return document
    if isinstance(document.get("reactions"), list):
        return {
            "schema_version": SCHEMA_VERSION,
            "workflow": WORKFLOW,
            "input_profile": "ord_dataset",
            "source": {
                "identifier": document.get("dataset_id") or path.name,
                "content_sha256": content_hash,
            },
            "options": {},
            "upstream_artifacts": [],
            "records": [
                {
                    "record_id": reaction.get("reaction_id")
                    or f"reaction-{index + 1}",
                    "ord_record": reaction,
                }
                for index, reaction in enumerate(document["reactions"])
                if isinstance(reaction, dict)
            ],
        }
    if "inputs" in document or "outcomes" in document:
        return {
            "schema_version": SCHEMA_VERSION,
            "workflow": WORKFLOW,
            "input_profile": "ord_reaction",
            "source": {
                "identifier": document.get("reaction_id") or path.name,
                "content_sha256": content_hash,
            },
            "options": {},
            "upstream_artifacts": [],
            "records": [
                {
                    "record_id": document.get("reaction_id") or "reaction-1",
                    "ord_record": document,
                }
            ],
        }
    raise InputFailure("无法识别 JSON 输入 profile")


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="输入 JSON/CSV/PB/PB.GZ")
    parser.add_argument("--output", required=True, help="输出 JSON")
    args = parser.parse_args(argv)
    try:
        toolkit = load_toolkit()
        request = load_request(Path(args.input), toolkit)
        result = process_request(request)
        output = Path(args.output)
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(
            json.dumps(result, ensure_ascii=False, indent=2) + "\n",
            encoding="utf-8",
        )
    except (OSError, InputFailure, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    print(
        f"完成 {result['input_summary']['output_records']} 条反应；"
        f"ready={result['input_summary']['disposition_counts']['ready_for_search']}，"
        f"review={result['input_summary']['disposition_counts']['review_required']}，"
        f"rejected={result['input_summary']['disposition_counts']['rejected']}。"
    )
    fatal_codes = {
        "E-INPUT-SCHEMA-001",
        "E-INPUT-HASH-001",
        "E-RESOURCE-LIMIT-001",
    }
    return (
        1
        if any(item.get("code") in fatal_codes for item in result["errors"])
        else 0
    )


if __name__ == "__main__":
    raise SystemExit(main())
