#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import html
import json
import math
import re
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


SCHEMA_VERSION = "1.0.0"
WORKFLOW = "predict-complex-structures"
REQUIRED_JSON = {
    "input_manifest.json": "input_fingerprint",
    "evidence.json": "evidence_fingerprint",
    "ligand_handoff.json": "ligand_handoff_fingerprint",
    "run_manifest.json": "run_fingerprint",
}
RUN_STATUSES = {"completed", "partial", "not_run", "error"}
REVIEW_DISPOSITIONS = {"ready_for_review", "review_required", "rejected"}
OFFICIAL_ENTITY_KEYS = {
    "proteinChain",
    "dnaSequence",
    "rnaSequence",
    "ligand",
    "ion",
}
TEXT_SAFETY_SUFFIXES = {
    ".a3m",
    ".cif",
    ".fa",
    ".fasta",
    ".hhr",
    ".html",
    ".json",
    ".mmcif",
    ".mol",
    ".mol2",
    ".pdb",
    ".sdf",
    ".txt",
}
STRUCTURE_NAME_RE = re.compile(
    r"^(?P<job>.+)_sample_(?P<sample>\d+)\.(?:cif|mmcif)$"
)
CONFIDENCE_NAME_RE = re.compile(
    r"^(?P<job>.+)_summary_confidence_sample_(?P<sample>\d+)\.json$"
)
FULL_CONFIDENCE_NAME_RE = re.compile(
    r"^(?P<job>.+)_full_data_sample_(?P<sample>\d+)\.json$"
)
SEED_DIR_RE = re.compile(r"^seed_(?P<seed>\d+)$")
SECRET_RE = re.compile(
    r"(?i)(?:[\"']?(?:api[_ -]?key|authorization|cookie|"
    r"x-agent-plan-key)[\"']?\s*[:=]\s*[\"']?[^\"'\s,}]{4,}|"
    r"bearer\s+[A-Za-z0-9._~-]{16,}|"
    r"ark-[A-Za-z0-9-]{30,}|"
    r"sk-[A-Za-z0-9_-]{20,}|"
    r"AKIA[0-9A-Z]{16}|"
    r"github_pat_[A-Za-z0-9_]{20,}|"
    r"gh[pousr]_[A-Za-z0-9]{20,}|"
    r"eyJ[A-Za-z0-9_-]{20,}\.[A-Za-z0-9_-]{20,}\.[A-Za-z0-9_-]{20,})"
)
SENSITIVE_PATH_RE = re.compile(
    r"(?:/"
    + "Users/"
    + r"[^/\s]+|/home/[^/\s]+|[A-Za-z]:\\"
    + "Users"
    + r"\\[^\\\s]+|/root(?:/[^\s]*)?|/mlx_devbox/users/[^/\s]+|"
    r"/private/(?:tmp|var)/[^\s]+)"
)
TEMPORAL_KEYS = {
    "generated_at_utc",
    "started_at_utc",
    "finished_at_utc",
    "retrieved_at_utc",
    "requested_at_utc",
    "confirmed_at_utc",
    "elapsed_seconds",
}
FINGERPRINT_KEYS = {
    "input_fingerprint",
    "evidence_fingerprint",
    "ligand_handoff_fingerprint",
    "run_fingerprint",
    "result_fingerprint",
}
JSON_LIMIT_BYTES = 64 * 1024 * 1024
STRUCTURE_LIMIT_BYTES = 512 * 1024 * 1024
GENERATED_REPORT_NAMES = {"validation_report.json", "report.html"}
PROTEIN_RESIDUE_NAMES = {
    "A": "ALA",
    "C": "CYS",
    "D": "ASP",
    "E": "GLU",
    "F": "PHE",
    "G": "GLY",
    "H": "HIS",
    "I": "ILE",
    "K": "LYS",
    "L": "LEU",
    "M": "MET",
    "N": "ASN",
    "P": "PRO",
    "Q": "GLN",
    "R": "ARG",
    "S": "SER",
    "T": "THR",
    "V": "VAL",
    "W": "TRP",
    "Y": "TYR",
    "X": "UNK",
}
DNA_RESIDUE_NAMES = {
    "A": "DA",
    "T": "DT",
    "G": "DG",
    "C": "DC",
    "N": "DN",
}
RNA_RESIDUE_NAMES = {
    "A": "A",
    "U": "U",
    "G": "G",
    "C": "C",
    "N": "N",
}


class ValidationDependencyFailure(RuntimeError):
    pass


def now_utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def canonical_json(value: Any) -> str:
    return json.dumps(
        value,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    )


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _without_nondeterministic_fields(value: Any) -> Any:
    if isinstance(value, dict):
        return {
            key: _without_nondeterministic_fields(item)
            for key, item in value.items()
            if key not in TEMPORAL_KEYS and key not in FINGERPRINT_KEYS
        }
    if isinstance(value, list):
        return [_without_nondeterministic_fields(item) for item in value]
    return value


def document_fingerprint(document: dict[str, Any]) -> str:
    return sha256_bytes(
        canonical_json(_without_nondeterministic_fields(document)).encode("utf-8")
    )


def read_json(path: Path, errors: list[str], label: str) -> Any | None:
    if path.is_symlink() or not path.is_file():
        errors.append(f"{label}: 文件缺失或不是普通文件")
        return None
    if path.stat().st_size > JSON_LIMIT_BYTES:
        errors.append(f"{label}: 超过 {JSON_LIMIT_BYTES} bytes")
        return None
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        errors.append(f"{label}: 无法读取有效 UTF-8 JSON: {exc}")
        return None


def write_json(path: Path, document: dict[str, Any]) -> None:
    write_text_safely(
        path,
        json.dumps(document, ensure_ascii=False, indent=2) + "\n",
    )


def write_text_safely(path: Path, text: str) -> None:
    if path.is_symlink():
        raise OSError(f"{path.name} 不允许 symlink")
    temporary = path.with_name(f".{path.name}.tmp")
    temporary.write_text(text, encoding="utf-8")
    temporary.replace(path)


def validate_fingerprint(
    document: Any,
    key: str,
    label: str,
    errors: list[str],
) -> None:
    if not isinstance(document, dict):
        return
    expected = document.get(key)
    if not isinstance(expected, str):
        errors.append(f"{label}: 缺少 {key}")
        return
    actual = document_fingerprint(document)
    if expected != actual:
        errors.append(f"{label}: {key} 与内容不一致")


def recursive_nonfinite(value: Any, path: str = "$") -> list[str]:
    findings: list[str] = []
    if isinstance(value, float) and not math.isfinite(value):
        findings.append(path)
    elif isinstance(value, dict):
        for key, item in value.items():
            findings.extend(recursive_nonfinite(item, f"{path}.{key}"))
    elif isinstance(value, list):
        for index, item in enumerate(value):
            findings.extend(recursive_nonfinite(item, f"{path}[{index}]"))
    return findings


def resolve_output_path(
    result_dir: Path,
    raw: Any,
    label: str,
    errors: list[str],
) -> Path | None:
    if not isinstance(raw, str) or not raw:
        errors.append(f"{label}: 路径必须是非空字符串")
        return None
    relative = Path(raw)
    if relative.is_absolute() or ".." in relative.parts:
        errors.append(f"{label}: 必须是结果目录内的相对路径")
        return None
    unresolved_candidate = result_dir / relative
    if unresolved_candidate.is_symlink():
        errors.append(f"{label}: 不允许 symlink")
        return None
    candidate = unresolved_candidate.resolve()
    try:
        candidate.relative_to(result_dir.resolve())
    except ValueError:
        errors.append(f"{label}: 路径逃逸结果目录")
        return None
    return candidate


def expected_names_from_payload(
    entity_key: str,
    payload: dict[str, Any],
) -> list[str] | None:
    sequence = payload.get("sequence")
    residue_map = {
        "proteinChain": PROTEIN_RESIDUE_NAMES,
        "dnaSequence": DNA_RESIDUE_NAMES,
        "rnaSequence": RNA_RESIDUE_NAMES,
    }.get(entity_key)
    if residue_map is not None and isinstance(sequence, str):
        names = [residue_map[residue] for residue in sequence]
        type_key, position_key = (
            ("ptmType", "ptmPosition")
            if entity_key == "proteinChain"
            else ("modificationType", "basePosition")
        )
        for modification in payload.get("modifications") or []:
            if not isinstance(modification, dict):
                continue
            position = modification.get(position_key)
            code = modification.get(type_key)
            if (
                isinstance(position, int)
                and 1 <= position <= len(names)
                and isinstance(code, str)
            ):
                names[position - 1] = code.removeprefix("CCD_")
        return ["MET" if name == "MSE" else name for name in names]
    if entity_key == "ion" and isinstance(payload.get("ion"), str):
        return [payload["ion"]]
    if entity_key == "ligand" and isinstance(payload.get("ligand"), str):
        ligand = payload["ligand"]
        if ligand.startswith("CCD_"):
            return ligand.removeprefix("CCD_").split("_")
    return None


def prediction_output_identity(
    relative_path: str,
    kind: str,
) -> dict[str, Any] | None:
    path = Path(relative_path)
    seeds = [
        int(match.group("seed"))
        for part in path.parts
        if (match := SEED_DIR_RE.fullmatch(part))
    ]
    if len(seeds) != 1:
        return None
    pattern = {
        "structure": STRUCTURE_NAME_RE,
        "confidence": CONFIDENCE_NAME_RE,
        "full_confidence": FULL_CONFIDENCE_NAME_RE,
    }.get(kind)
    if pattern is None:
        return None
    match = pattern.fullmatch(path.name)
    if match is None:
        return None
    job = match.group("job")
    sample = int(match.group("sample"))
    return {
        "job": job,
        "seed": seeds[0],
        "sample": sample,
        "key": f"{job}|{seeds[0]}|{sample}",
    }


def validate_protenix_input(
    document: Any,
    input_manifest: dict[str, Any],
    run_manifest: dict[str, Any],
    errors: list[str],
) -> dict[str, Any]:
    summary = {
        "task_count": 0,
        "entity_count": 0,
        "chain_ids": [],
        "chain_expectations": {},
    }
    if not isinstance(document, list):
        errors.append("protenix_input.json: 顶层必须是 array")
        return summary
    if not run_manifest.get("inference_allowed"):
        if document:
            errors.append("protenix_input.json: 上游阻断时必须为空数组")
        return summary
    if not document:
        errors.append("protenix_input.json: 可执行任务不得为空")
        return summary
    chain_ids: list[str] = []
    for task_index, task in enumerate(document):
        if not isinstance(task, dict):
            errors.append(f"protenix_input[{task_index}]: 必须是 object")
            continue
        unknown = set(task) - {"name", "sequences", "covalent_bonds", "constraint"}
        if unknown:
            errors.append(
                f"protenix_input[{task_index}]: 非官方顶层字段 {sorted(unknown)}"
            )
        if task.get("name") != input_manifest.get("job_id"):
            errors.append(f"protenix_input[{task_index}].name 与 job_id 不一致")
        sequences = task.get("sequences")
        if not isinstance(sequences, list) or not sequences:
            errors.append(f"protenix_input[{task_index}].sequences 不能为空")
            continue
        for entity_index, entity in enumerate(sequences):
            if not isinstance(entity, dict) or len(entity) != 1:
                errors.append(
                    f"protenix_input[{task_index}].sequences[{entity_index}] "
                    "必须恰好包含一个实体类型"
                )
                continue
            entity_key = next(iter(entity))
            if entity_key not in OFFICIAL_ENTITY_KEYS:
                errors.append(
                    f"protenix_input[{task_index}].sequences[{entity_index}] "
                    f"实体类型不受支持: {entity_key}"
                )
                continue
            payload = entity[entity_key]
            if not isinstance(payload, dict):
                errors.append(
                    f"protenix_input[{task_index}].sequences[{entity_index}] "
                    "实体值必须是 object"
                )
                continue
            ids = payload.get("id")
            if not isinstance(ids, list) or not ids:
                errors.append(
                    f"protenix_input[{task_index}].sequences[{entity_index}]."
                    f"{entity_key}.id 必须是非空 array"
                )
            else:
                chain_ids.extend(str(item) for item in ids)
                count = payload.get("count")
                if (
                    isinstance(count, bool)
                    or not isinstance(count, int)
                    or count <= 0
                    or len(ids) != count
                ):
                    errors.append(
                        f"protenix_input[{task_index}].sequences[{entity_index}]."
                        f"{entity_key}: count 必须与 id 数量一致"
                    )
                sequence = payload.get("sequence")
                expected_residue_count = (
                    len(sequence)
                    if entity_key
                    in {"proteinChain", "dnaSequence", "rnaSequence"}
                    and isinstance(sequence, str)
                    else None
                )
                expected_residue_names = expected_names_from_payload(
                    entity_key,
                    payload,
                )
                if expected_residue_count is None and expected_residue_names:
                    expected_residue_count = len(expected_residue_names)
                for chain_id in ids:
                    summary["chain_expectations"][str(chain_id)] = {
                        "entity_type": entity_key,
                        "residue_count": expected_residue_count,
                        "residue_names": expected_residue_names,
                    }
            summary["entity_count"] += 1
    expected_ids = set(input_manifest.get("entity_chain_ids") or [])
    if set(chain_ids) != expected_ids:
        errors.append(
            "protenix_input.json: 实体链 ID 与 input_manifest 不一致"
        )
    summary["task_count"] = len(document)
    summary["chain_ids"] = sorted(chain_ids)
    return summary


def load_gemmi() -> Any:
    try:
        import gemmi
    except ImportError as exc:
        raise ValidationDependencyFailure(
            "缺少 gemmi==0.6.7，无法独立解析 CIF/mmCIF。"
        ) from exc
    version = getattr(gemmi, "__version__", None)
    if version != "0.6.7":
        raise ValidationDependencyFailure(
            f"需要 gemmi==0.6.7，当前为 {version}。"
        )
    return gemmi


def parse_structure(
    path: Path,
    gemmi: Any,
) -> dict[str, Any]:
    try:
        structure = gemmi.read_structure(str(path))
    except Exception as exc:
        raise ValueError(f"mmCIF 无法解析: {exc}") from exc
    model_count = len(structure)
    if model_count < 1:
        raise ValueError("mmCIF 不包含 model")
    chain_profiles: dict[str, dict[str, Any]] = {}
    residue_count = 0
    atom_count = 0
    for model in structure:
        for chain in model:
            chain_name = str(chain.name)
            profile = chain_profiles.setdefault(
                chain_name,
                {
                    "residue_names": [],
                    "heavy_atom_count": 0,
                    "heavy_element_counts": Counter(),
                },
            )
            for residue in chain:
                residue_count += 1
                profile["residue_names"].append(str(residue.name))
                for atom in residue:
                    atom_count += 1
                    element = str(atom.element.name).strip()
                    if element and element.upper() != "H":
                        profile["heavy_atom_count"] += 1
                        profile["heavy_element_counts"][element] += 1
    if residue_count < 1 or atom_count < 1:
        raise ValueError("mmCIF 不包含 residue/atom")
    normalized_profiles = {}
    for chain_id, profile in chain_profiles.items():
        normalized_profiles[chain_id] = {
            "residue_names": profile["residue_names"],
            "heavy_atom_count": profile["heavy_atom_count"],
            "heavy_element_counts": dict(
                sorted(profile["heavy_element_counts"].items())
            ),
        }
    return {
        "model_count": model_count,
        "chain_ids": sorted(chain_profiles),
        "residue_count": residue_count,
        "atom_count": atom_count,
        "chain_profiles": normalized_profiles,
    }


def validate_ligand_structure(
    structure_profile: dict[str, Any],
    ligand_handoff: dict[str, Any],
    path: str,
    errors: list[str],
    warnings: list[str],
) -> dict[str, Any]:
    result = {
        "file": path,
        "status": "not_applicable",
        "chain_id": ligand_handoff.get("chain_id"),
    }
    if ligand_handoff.get("ligand_status") == "not_applicable":
        return result
    chain_id = ligand_handoff.get("chain_id")
    chain_profile = structure_profile["chain_profiles"].get(chain_id)
    if chain_profile is None:
        errors.append(f"{path}: 缺少配体链 {chain_id}")
        result["status"] = "error"
        return result
    expected = ligand_handoff.get("expected_output_profile")
    if not isinstance(expected, dict):
        warnings.append(f"{path}: 没有可独立核对的配体原子/CCD 口径")
        result["status"] = "not_assessed"
        return result
    if expected.get("ccd_id"):
        expected_names = str(expected["ccd_id"]).removeprefix("CCD_").split("_")
        if chain_profile["residue_names"] != expected_names:
            errors.append(
                f"{path}: 配体链 residue 与 CCD {expected['ccd_id']} 不一致"
            )
            result["status"] = "error"
            return result
        result["status"] = "matched_ccd"
        return result
    expected_count = expected.get("heavy_atom_count")
    expected_elements = expected.get("heavy_element_counts")
    if chain_profile["heavy_atom_count"] != expected_count:
        errors.append(
            f"{path}: 配体重原子数不一致，expected={expected_count} "
            f"actual={chain_profile['heavy_atom_count']}"
        )
    if chain_profile["heavy_element_counts"] != expected_elements:
        errors.append(
            f"{path}: 配体重元素组成与已审查结构不一致"
        )
    result.update(
        {
            "expected_heavy_atom_count": expected_count,
            "actual_heavy_atom_count": chain_profile["heavy_atom_count"],
            "expected_heavy_element_counts": expected_elements,
            "actual_heavy_element_counts": chain_profile[
                "heavy_element_counts"
            ],
            "status": (
                "matched"
                if chain_profile["heavy_atom_count"] == expected_count
                and chain_profile["heavy_element_counts"] == expected_elements
                else "error"
            ),
        }
    )
    return result


def validate_structure_files(
    result_dir: Path,
    run_manifest: dict[str, Any],
    input_manifest: dict[str, Any],
    ligand_handoff: dict[str, Any],
    chain_expectations: dict[str, dict[str, Any]],
    errors: list[str],
    warnings: list[str],
) -> list[dict[str, Any]]:
    outputs = run_manifest.get("outputs")
    raw_paths = (
        outputs.get("structure_files")
        if isinstance(outputs, dict)
        else None
    )
    if not isinstance(raw_paths, list):
        errors.append("run_manifest.outputs.structure_files 必须是 array")
        return []
    if run_manifest.get("status") in {"completed", "partial"} and not raw_paths:
        errors.append("Protenix 输出缺少 CIF/mmCIF")
        return []
    if not raw_paths:
        return []
    gemmi = load_gemmi()
    expected_chain_ids = set(input_manifest.get("entity_chain_ids") or [])
    checks: list[dict[str, Any]] = []
    for index, raw_path in enumerate(raw_paths):
        output_identity = prediction_output_identity(raw_path, "structure")
        if output_identity is None:
            errors.append(
                f"{raw_path}: 无法从结构路径解析唯一 job/seed/sample"
            )
        path = resolve_output_path(
            result_dir,
            raw_path,
            f"structure_files[{index}]",
            errors,
        )
        if path is None:
            continue
        if not path.is_file() or path.suffix.lower() not in {".cif", ".mmcif"}:
            errors.append(f"{raw_path}: CIF/mmCIF 文件缺失或扩展名错误")
            continue
        if path.stat().st_size > STRUCTURE_LIMIT_BYTES:
            errors.append(
                f"{raw_path}: 超过 {STRUCTURE_LIMIT_BYTES} bytes 结构文件上限"
            )
            continue
        try:
            profile = parse_structure(path, gemmi)
        except ValueError as exc:
            errors.append(f"{raw_path}: {exc}")
            checks.append({"file": raw_path, "status": "error"})
            continue
        actual_chain_ids = set(profile["chain_ids"])
        missing = sorted(expected_chain_ids - actual_chain_ids)
        unexpected = sorted(actual_chain_ids - expected_chain_ids)
        unexpected_model_count = profile["model_count"] != 1
        if missing:
            errors.append(f"{raw_path}: 缺少期望链 {missing}")
        if unexpected:
            errors.append(f"{raw_path}: 包含未请求链 {unexpected}")
        if unexpected_model_count:
            errors.append(
                f"{raw_path}: 每个预测样本必须恰好包含 1 个 model，"
                f"actual={profile['model_count']}"
            )
        residue_mismatches = []
        sequence_mismatches = []
        for chain_id, expectation in chain_expectations.items():
            expected_count = expectation.get("residue_count")
            if expected_count is None or chain_id not in profile["chain_profiles"]:
                continue
            actual_count = len(
                profile["chain_profiles"][chain_id]["residue_names"]
            )
            if actual_count != expected_count:
                residue_mismatches.append(
                    {
                        "chain_id": chain_id,
                        "expected": expected_count,
                        "actual": actual_count,
                    }
                )
                errors.append(
                    f"{raw_path}: 链 {chain_id} 残基数不一致，"
                    f"expected={expected_count} actual={actual_count}"
                )
            expected_names = expectation.get("residue_names")
            actual_names = profile["chain_profiles"][chain_id][
                "residue_names"
            ]
            if (
                isinstance(expected_names, list)
                and len(actual_names) == len(expected_names)
                and actual_names != expected_names
            ):
                mismatch_positions = [
                    index + 1
                    for index, (expected, actual) in enumerate(
                        zip(expected_names, actual_names)
                    )
                    if expected != actual
                ]
                sequence_mismatches.append(
                    {
                        "chain_id": chain_id,
                        "positions": mismatch_positions,
                    }
                )
                errors.append(
                    f"{raw_path}: 链 {chain_id} 残基序列与输入不一致，"
                    f"positions={mismatch_positions[:20]}"
                )
        ligand_check = validate_ligand_structure(
            profile,
            ligand_handoff,
            raw_path,
            errors,
            warnings,
        )
        checks.append(
            {
                "file": raw_path,
                "sha256": sha256_file(path),
                "status": (
                    "error"
                    if output_identity is None
                    or missing
                    or unexpected
                    or unexpected_model_count
                    or residue_mismatches
                    or sequence_mismatches
                    or ligand_check["status"] == "error"
                    else "valid"
                ),
                "output_identity": output_identity,
                "profile": profile,
                "ligand_check": ligand_check,
                "residue_mismatches": residue_mismatches,
                "sequence_mismatches": sequence_mismatches,
            }
        )
    return checks


def validate_auxiliary_structure_files(
    result_dir: Path,
    run_manifest: dict[str, Any],
    errors: list[str],
) -> list[dict[str, Any]]:
    outputs = run_manifest.get("outputs")
    raw_paths = (
        outputs.get("auxiliary_structure_files")
        if isinstance(outputs, dict)
        else None
    )
    if not isinstance(raw_paths, list):
        errors.append("run_manifest.outputs.auxiliary_structure_files 必须是 array")
        return []
    if not raw_paths:
        return []
    gemmi = load_gemmi()
    checks: list[dict[str, Any]] = []
    for index, raw_path in enumerate(raw_paths):
        path = resolve_output_path(
            result_dir,
            raw_path,
            f"auxiliary_structure_files[{index}]",
            errors,
        )
        if path is None:
            continue
        if (
            not path.is_file()
            or path.suffix.lower() not in {".cif", ".mmcif"}
            or not path.stem.endswith("_wounresol")
        ):
            errors.append(f"{raw_path}: 辅助结构文件名或扩展名不合法")
            continue
        if path.stat().st_size > STRUCTURE_LIMIT_BYTES:
            errors.append(f"{raw_path}: 超过结构文件上限")
            continue
        try:
            profile = parse_structure(path, gemmi)
        except ValueError as exc:
            errors.append(f"{raw_path}: {exc}")
            checks.append({"file": raw_path, "status": "error"})
            continue
        checks.append(
            {
                "file": raw_path,
                "status": "valid",
                "sha256": sha256_file(path),
                "profile": profile,
            }
        )
    return checks


def validate_confidence_document(
    document: Any,
    path: str,
    complex_task: bool,
    chain_count: int,
    errors: list[str],
) -> dict[str, Any]:
    result = {"file": path, "status": "valid", "keys": []}
    if not isinstance(document, dict):
        errors.append(f"{path}: summary confidence 顶层必须是 object")
        result["status"] = "error"
        return result
    initial_error_count = len(errors)
    required = {"ranking_score", "plddt", "ptm"}
    missing = sorted(required - set(document))
    if missing:
        errors.append(f"{path}: 缺少 confidence 字段 {missing}")
    if complex_task and "iptm" not in document:
        errors.append(f"{path}: 复合物任务缺少 iptm")
    if complex_task and not (
        "chain_pair_iptm" in document
        or "chain_pair_iptm_global" in document
    ):
        errors.append(f"{path}: 复合物任务缺少 chain-pair ipTM")
    nonfinite = recursive_nonfinite(document)
    if nonfinite:
        errors.append(f"{path}: 含非有限数值 {nonfinite[:10]}")

    def validate_scalar(key: str, minimum: float | None, maximum: float | None):
        value = document.get(key)
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            errors.append(f"{path}: {key} 必须是有限数值")
            return
        numeric = float(value)
        if not math.isfinite(numeric):
            return
        if minimum is not None and numeric < minimum:
            errors.append(f"{path}: {key} 低于 {minimum}")
        if maximum is not None and numeric > maximum:
            errors.append(f"{path}: {key} 高于 {maximum}")

    validate_scalar("ranking_score", None, None)
    validate_scalar("plddt", 0.0, 100.0)
    validate_scalar("ptm", 0.0, 1.0)
    if complex_task:
        validate_scalar("iptm", 0.0, 1.0)
    if "gpde" in document:
        validate_scalar("gpde", 0.0, None)
    if "disorder" in document:
        validate_scalar("disorder", 0.0, 1.0)
    if "has_clash" in document and not isinstance(
        document["has_clash"], bool
    ):
        errors.append(f"{path}: has_clash 必须是 boolean")
    if "num_recycles" in document:
        num_recycles = document["num_recycles"]
        if (
            isinstance(num_recycles, bool)
            or not isinstance(num_recycles, int)
            or num_recycles < 0
        ):
            errors.append(f"{path}: num_recycles 必须是非负整数")

    def validate_vector(key: str, minimum: float, maximum: float) -> None:
        if key not in document:
            return
        vector = document[key]
        if not isinstance(vector, list) or len(vector) != chain_count:
            errors.append(f"{path}: {key} 必须是长度 {chain_count} 的数组")
            return
        for index, value in enumerate(vector):
            if (
                isinstance(value, bool)
                or not isinstance(value, (int, float))
                or not math.isfinite(float(value))
                or not minimum <= float(value) <= maximum
            ):
                errors.append(
                    f"{path}: {key}[{index}] 必须是 {minimum}–{maximum} 有限数值"
                )

    validate_vector("chain_ptm", 0.0, 1.0)
    validate_vector("chain_iptm", 0.0, 1.0)
    validate_vector("chain_plddt", 0.0, 100.0)

    matrix_ranges = {
        "chain_pair_iptm": (0.0, 1.0),
        "chain_pair_iptm_global": (0.0, 1.0),
        "chain_pair_plddt": (0.0, 100.0),
    }
    for matrix_key, (minimum, maximum) in matrix_ranges.items():
        if matrix_key not in document:
            continue
        matrix = document[matrix_key]
        if (
            not isinstance(matrix, list)
            or len(matrix) != chain_count
            or any(
                not isinstance(row, list) or len(row) != chain_count
                for row in matrix
            )
        ):
            errors.append(
                f"{path}: {matrix_key} 必须是 {chain_count}x{chain_count} 矩阵"
            )
            continue
        for row_index, row in enumerate(matrix):
            for column_index, value in enumerate(row):
                if (
                    isinstance(value, bool)
                    or not isinstance(value, (int, float))
                    or not math.isfinite(float(value))
                    or not minimum <= float(value) <= maximum
                ):
                    errors.append(
                        f"{path}: {matrix_key}[{row_index}]"
                        f"[{column_index}] 必须是 {minimum}–{maximum} 有限数值"
                    )

    result["keys"] = sorted(document)
    if len(errors) > initial_error_count:
        result["status"] = "error"
    return result


def validate_confidence_files(
    result_dir: Path,
    run_manifest: dict[str, Any],
    input_manifest: dict[str, Any],
    errors: list[str],
) -> list[dict[str, Any]]:
    outputs = run_manifest.get("outputs")
    raw_paths = (
        outputs.get("summary_confidence_files")
        if isinstance(outputs, dict)
        else None
    )
    if not isinstance(raw_paths, list):
        errors.append(
            "run_manifest.outputs.summary_confidence_files 必须是 array"
        )
        return []
    if run_manifest.get("status") in {"completed", "partial"} and not raw_paths:
        errors.append("Protenix 输出缺少 summary confidence JSON")
        return []
    complex_task = len(input_manifest.get("entity_chain_ids") or []) > 1
    chain_count = len(input_manifest.get("entity_chain_ids") or [])
    checks: list[dict[str, Any]] = []
    for index, raw_path in enumerate(raw_paths):
        output_identity = prediction_output_identity(raw_path, "confidence")
        if output_identity is None:
            errors.append(
                f"{raw_path}: 无法从 confidence 路径解析唯一 job/seed/sample"
            )
        path = resolve_output_path(
            result_dir,
            raw_path,
            f"summary_confidence_files[{index}]",
            errors,
        )
        if path is None:
            continue
        document = read_json(path, errors, raw_path)
        check = validate_confidence_document(
            document,
            raw_path,
            complex_task,
            chain_count,
            errors,
        )
        if path.is_file():
            check["sha256"] = sha256_file(path)
        check["output_identity"] = output_identity
        if output_identity is None:
            check["status"] = "error"
        checks.append(check)
    return checks


def validate_full_confidence_files(
    result_dir: Path,
    run_manifest: dict[str, Any],
    errors: list[str],
) -> list[dict[str, Any]]:
    outputs = run_manifest.get("outputs")
    raw_paths = (
        outputs.get("full_confidence_files")
        if isinstance(outputs, dict)
        else None
    )
    if not isinstance(raw_paths, list):
        errors.append("run_manifest.outputs.full_confidence_files 必须是 array")
        return []
    checks: list[dict[str, Any]] = []
    for index, raw_path in enumerate(raw_paths):
        identity = prediction_output_identity(raw_path, "full_confidence")
        if identity is None:
            errors.append(
                f"{raw_path}: 无法从 full confidence 路径解析唯一 job/seed/sample"
            )
        path = resolve_output_path(
            result_dir,
            raw_path,
            f"full_confidence_files[{index}]",
            errors,
        )
        if path is None:
            continue
        document = read_json(path, errors, raw_path)
        nonfinite = recursive_nonfinite(document)
        if nonfinite:
            errors.append(f"{raw_path}: 含非有限数值 {nonfinite[:10]}")
        checks.append(
            {
                "file": raw_path,
                "status": (
                    "valid"
                    if identity is not None
                    and document is not None
                    and not nonfinite
                    else "error"
                ),
                "sha256": sha256_file(path) if path.is_file() else None,
                "output_identity": identity,
            }
        )
    return checks


def validate_output_pairing(
    structure_checks: list[dict[str, Any]],
    confidence_checks: list[dict[str, Any]],
    run_manifest: dict[str, Any],
    errors: list[str],
) -> dict[str, Any]:
    structure_keys = [
        item["output_identity"]["key"]
        for item in structure_checks
        if isinstance(item.get("output_identity"), dict)
    ]
    confidence_keys = [
        item["output_identity"]["key"]
        for item in confidence_checks
        if isinstance(item.get("output_identity"), dict)
    ]
    if len(structure_keys) != len(set(structure_keys)):
        errors.append("Protenix 结构输出包含重复 job/seed/sample key")
    if len(confidence_keys) != len(set(confidence_keys)):
        errors.append("Protenix confidence 输出包含重复 job/seed/sample key")
    structure_set = set(structure_keys)
    confidence_set = set(confidence_keys)
    missing_confidence = sorted(structure_set - confidence_set)
    missing_structure = sorted(confidence_set - structure_set)
    if missing_confidence or missing_structure:
        errors.append(
            "Protenix 结构与 summary confidence 的 job/seed/sample 不一致"
        )
    model = run_manifest.get("model")
    requested_seeds = (
        set(model.get("seeds") or []) if isinstance(model, dict) else set()
    )
    requested_samples = (
        model.get("samples") if isinstance(model, dict) else None
    )
    job_id = run_manifest.get("job_id")
    expected_keys = (
        {
            f"{job_id}|{seed}|{sample}"
            for seed in requested_seeds
            for sample in range(requested_samples)
        }
        if isinstance(requested_samples, int)
        and requested_samples > 0
        and isinstance(job_id, str)
        else set()
    )
    observed_seeds = {
        item["output_identity"]["seed"]
        for item in structure_checks + confidence_checks
        if isinstance(item.get("output_identity"), dict)
    }
    unexpected_seeds = sorted(observed_seeds - requested_seeds)
    if unexpected_seeds:
        errors.append(f"Protenix 输出包含未请求的 seed：{unexpected_seeds}")
    missing_expected = sorted(expected_keys - structure_set)
    unexpected_output_keys = sorted(structure_set - expected_keys)
    if missing_expected:
        errors.append(
            f"Protenix 输出缺少请求的 job/seed/sample：{missing_expected}"
        )
    if unexpected_output_keys:
        errors.append(
            f"Protenix 输出包含未请求的 job/seed/sample："
            f"{unexpected_output_keys}"
        )
    return {
        "status": (
            "matched"
            if not missing_confidence
            and not missing_structure
            and not unexpected_seeds
            and not missing_expected
            and not unexpected_output_keys
            and len(structure_keys) == len(structure_set)
            and len(confidence_keys) == len(confidence_set)
            else "error"
        ),
        "structure_keys": sorted(structure_set),
        "confidence_keys": sorted(confidence_set),
        "missing_confidence": missing_confidence,
        "missing_structure": missing_structure,
        "requested_seeds": sorted(requested_seeds),
        "requested_samples_per_seed": requested_samples,
        "expected_keys": sorted(expected_keys),
        "observed_seeds": sorted(observed_seeds),
        "unexpected_seeds": unexpected_seeds,
        "missing_expected": missing_expected,
        "unexpected_output_keys": unexpected_output_keys,
    }


def scan_result_safety(result_dir: Path, errors: list[str]) -> list[str]:
    scanned: list[str] = []
    for path in sorted(result_dir.rglob("*")):
        if path.is_symlink():
            errors.append(
                f"{path.relative_to(result_dir)}: 结果目录不允许 symlink"
            )
            continue
        if not path.is_file():
            continue
        if path.name in GENERATED_REPORT_NAMES:
            continue
        if path.suffix.lower() not in TEXT_SAFETY_SUFFIXES:
            continue
        try:
            handle = path.open(encoding="utf-8", errors="replace")
        except OSError as exc:
            errors.append(
                f"{path.relative_to(result_dir)}: 安全扫描无法读取文件: {exc}"
            )
            continue
        relative = path.relative_to(result_dir).as_posix()
        scanned.append(relative)
        with handle:
            for line_number, line in enumerate(handle, start=1):
                if SECRET_RE.search(line):
                    errors.append(
                        f"{relative}:{line_number}: possible credential detected"
                    )
                if SENSITIVE_PATH_RE.search(line):
                    errors.append(
                        f"{relative}:{line_number}: host-sensitive path detected"
                    )
    return scanned


def discover_declared_outputs(result_dir: Path) -> dict[str, list[str]]:
    output_root = result_dir / "protenix_output"
    result = {
        "structure_files": [],
        "auxiliary_structure_files": [],
        "summary_confidence_files": [],
        "full_confidence_files": [],
    }
    if not output_root.is_dir():
        return result
    for path in sorted(output_root.rglob("*")):
        if path.is_symlink() or not path.is_file():
            continue
        relative = path.relative_to(result_dir).as_posix()
        if path.suffix.lower() in {".cif", ".mmcif"}:
            if path.stem.endswith("_wounresol"):
                result["auxiliary_structure_files"].append(relative)
            else:
                result["structure_files"].append(relative)
        elif path.suffix.lower() == ".json" and (
            "_summary_confidence_sample_" in path.name
        ):
            result["summary_confidence_files"].append(relative)
        elif path.suffix.lower() == ".json" and (
            "_full_data_sample_" in path.name
        ):
            result["full_confidence_files"].append(relative)
    return result


def validate_output_inventory(
    result_dir: Path,
    run_manifest: dict[str, Any],
    errors: list[str],
) -> dict[str, Any]:
    actual = discover_declared_outputs(result_dir)
    declared_root = run_manifest.get("outputs")
    declared = declared_root if isinstance(declared_root, dict) else {}
    mismatches: dict[str, dict[str, list[str]]] = {}
    for key, actual_paths in actual.items():
        declared_paths = declared.get(key)
        if not isinstance(declared_paths, list):
            continue
        missing_from_manifest = sorted(set(actual_paths) - set(declared_paths))
        missing_on_disk = sorted(set(declared_paths) - set(actual_paths))
        if missing_from_manifest or missing_on_disk:
            mismatches[key] = {
                "missing_from_manifest": missing_from_manifest,
                "missing_on_disk": missing_on_disk,
            }
            errors.append(
                f"run_manifest.outputs.{key} 与磁盘实际输出不一致"
            )
    return {"actual": actual, "mismatches": mismatches}


def render_report(document: dict[str, Any]) -> str:
    status = html.escape(str(document.get("status")))
    disposition = html.escape(str(document.get("review_disposition")))
    errors = document.get("errors") or []
    warnings = document.get("warnings") or []
    file_rows = []
    for item in document.get("file_checks", []):
        file_rows.append(
            "<tr><td>"
            + html.escape(str(item.get("file")))
            + "</td><td>"
            + html.escape(str(item.get("status")))
            + "</td></tr>"
        )
    issue_items = "".join(
        f"<li>{html.escape(str(item))}</li>" for item in errors + warnings
    )
    boundaries = "".join(
        f"<li>{html.escape(str(item))}</li>"
        for item in document.get("scientific_boundaries", [])
    )
    return f"""<!doctype html>
<html lang="zh-CN">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>Protenix 结构预测复核报告</title>
  <style>
    body {{ font-family: system-ui, sans-serif; max-width: 960px; margin: 40px auto;
      padding: 0 24px; color: #172033; background: #f7f8fb; line-height: 1.6; }}
    main {{ background: white; padding: 32px; border: 1px solid #e4e8f0;
      border-radius: 12px; }}
    h1, h2 {{ color: #0f172a; }}
    .status {{ display: grid; grid-template-columns: repeat(2, minmax(0, 1fr));
      gap: 16px; }}
    .card {{ background: #f0f2f7; padding: 16px; border-radius: 8px; }}
    table {{ width: 100%; border-collapse: collapse; }}
    th, td {{ text-align: left; border-bottom: 1px solid #e4e8f0; padding: 10px; }}
    .notice {{ border-left: 4px solid #d97706; padding-left: 16px; }}
  </style>
</head>
<body>
<main>
  <h1>生物分子复合物结构预测复核报告</h1>
  <p class="notice">该结果是模型生成的结构假设，不是实验结构、结合证明、
  活性结论或实验执行建议。</p>
  <section class="status">
    <div class="card"><strong>运行状态</strong><br>{status}</div>
    <div class="card"><strong>审查结论</strong><br>{disposition}</div>
  </section>
  <h2>文件检查</h2>
  <table><thead><tr><th>文件</th><th>状态</th></tr></thead>
  <tbody>{''.join(file_rows)}</tbody></table>
  <h2>错误与复核项</h2>
  <ul>{issue_items or '<li>无结构化错误；仍需专家复核。</li>'}</ul>
  <h2>科学边界</h2>
  <ul>{boundaries}</ul>
</main>
</body>
</html>
"""


def validate_result_dir(
    result_dir: Path,
    *,
    generated_at_utc: str | None = None,
    write_artifacts: bool = True,
) -> dict[str, Any]:
    if result_dir.is_symlink():
        raise OSError("result_dir 不允许 symlink")
    result_dir = result_dir.resolve()
    errors: list[str] = []
    warnings: list[str] = []
    notices: list[str] = []
    documents: dict[str, Any] = {}
    for name, fingerprint_key in REQUIRED_JSON.items():
        document = read_json(result_dir / name, errors, name)
        documents[name] = document
        validate_fingerprint(document, fingerprint_key, name, errors)
    protenix_input = read_json(
        result_dir / "protenix_input.json",
        errors,
        "protenix_input.json",
    )
    documents["protenix_input.json"] = protenix_input

    input_manifest = documents.get("input_manifest.json")
    evidence = documents.get("evidence.json")
    ligand_handoff = documents.get("ligand_handoff.json")
    run_manifest = documents.get("run_manifest.json")
    if not all(
        isinstance(item, dict)
        for item in (
            input_manifest,
            evidence,
            ligand_handoff,
            run_manifest,
        )
    ):
        input_manifest = input_manifest if isinstance(input_manifest, dict) else {}
        evidence = evidence if isinstance(evidence, dict) else {}
        ligand_handoff = (
            ligand_handoff if isinstance(ligand_handoff, dict) else {}
        )
        run_manifest = run_manifest if isinstance(run_manifest, dict) else {}

    if input_manifest.get("workflow") != WORKFLOW:
        errors.append("input_manifest.json: workflow 不正确")
    if run_manifest.get("workflow") != WORKFLOW:
        errors.append("run_manifest.json: workflow 不正确")
    status = run_manifest.get("status")
    if status not in RUN_STATUSES:
        errors.append(f"run_manifest.status 不受控: {status}")
    disposition = run_manifest.get("review_disposition")
    if disposition not in REVIEW_DISPOSITIONS:
        errors.append(
            f"run_manifest.review_disposition 不受控: {disposition}"
        )
    if (
        run_manifest.get("input_fingerprint")
        != input_manifest.get("input_fingerprint")
    ):
        errors.append("run_manifest.input_fingerprint 与 input_manifest 不一致")
    if run_manifest.get("job_id") != input_manifest.get("job_id"):
        errors.append("run_manifest.job_id 与 input_manifest 不一致")
    if run_manifest.get("model") != input_manifest.get("model_request"):
        errors.append("run_manifest.model 与 input_manifest 未绑定")
    actual_protenix_input_sha256 = sha256_bytes(
        canonical_json(protenix_input).encode("utf-8")
    )
    if (
        input_manifest.get("protenix_input_sha256")
        != actual_protenix_input_sha256
    ):
        errors.append("protenix_input.json 与 input_manifest 未绑定")
    if (
        input_manifest.get("bound_evidence_fingerprint")
        != evidence.get("evidence_fingerprint")
    ):
        errors.append("evidence.json 与 input_manifest 未绑定")
    if (
        input_manifest.get("bound_ligand_handoff_fingerprint")
        != ligand_handoff.get("ligand_handoff_fingerprint")
    ):
        errors.append("ligand_handoff.json 与 input_manifest 未绑定")
    if evidence.get("evidence_status") == "not_run":
        warnings.append("专业数据集证据未运行；不得伪造候选文献。")
    response = evidence.get("response")
    if isinstance(response, dict):
        for index, item in enumerate(response.get("items") or []):
            nested_doi = (
                item.get("extra_data", {}).get("doi")
                if isinstance(item, dict)
                and isinstance(item.get("extra_data"), dict)
                else ""
            )
            if (
                isinstance(item, dict)
                and not str(item.get("doi") or nested_doi or "").strip()
            ):
                warnings.append(f"evidence item {index} 缺少 DOI，原样保留。")
    if ligand_handoff.get("review_disposition") == "review_required":
        warnings.append("上游配体状态需要人工复核，风险已传播。")
    expected_scientific_review_reasons: list[str] = []
    if input_manifest.get("covalent_bonds"):
        expected_scientific_review_reasons.append(
            "COVALENT_BOND_NOT_INDEPENDENTLY_VERIFIED"
        )
    if input_manifest.get("constraint") is not None:
        expected_scientific_review_reasons.append(
            "SOFT_CONSTRAINT_NOT_INDEPENDENTLY_VERIFIED"
        )
    scientific_review_reasons = input_manifest.get(
        "scientific_review_reasons",
        [],
    )
    if not isinstance(scientific_review_reasons, list) or any(
        not isinstance(item, str) for item in scientific_review_reasons
    ):
        errors.append("input_manifest.scientific_review_reasons 必须是字符串数组")
        scientific_review_reasons = []
    if scientific_review_reasons != expected_scientific_review_reasons:
        errors.append(
            "input_manifest.scientific_review_reasons 与共价键/软约束输入不一致"
        )
    if run_manifest.get("scientific_review_reasons", []) != scientific_review_reasons:
        errors.append("run_manifest.scientific_review_reasons 与 input_manifest 不一致")
    if "COVALENT_BOND_NOT_INDEPENDENTLY_VERIFIED" in scientific_review_reasons:
        warnings.append("共价键意图尚未由独立 Validator 验证，必须人工复核。")
    if "SOFT_CONSTRAINT_NOT_INDEPENDENTLY_VERIFIED" in scientific_review_reasons:
        warnings.append("软约束满足情况尚未由独立 Validator 验证，必须人工复核。")

    input_check = validate_protenix_input(
        protenix_input,
        input_manifest,
        run_manifest,
        errors,
    )
    structure_checks: list[dict[str, Any]] = []
    auxiliary_structure_checks: list[dict[str, Any]] = []
    confidence_checks: list[dict[str, Any]] = []
    full_confidence_checks: list[dict[str, Any]] = []
    output_pairing = {
        "status": "not_run",
        "structure_keys": [],
        "confidence_keys": [],
    }
    output_inventory = validate_output_inventory(
        result_dir,
        run_manifest,
        errors,
    )
    if status in {"completed", "partial"}:
        try:
            structure_checks = validate_structure_files(
                result_dir,
                run_manifest,
                input_manifest,
                ligand_handoff,
                input_check.get("chain_expectations", {}),
                errors,
                warnings,
            )
        except ValidationDependencyFailure as exc:
            errors.append(str(exc))
        try:
            auxiliary_structure_checks = validate_auxiliary_structure_files(
                result_dir,
                run_manifest,
                errors,
            )
        except ValidationDependencyFailure as exc:
            errors.append(str(exc))
        confidence_checks = validate_confidence_files(
            result_dir,
            run_manifest,
            input_manifest,
            errors,
        )
        full_confidence_checks = validate_full_confidence_files(
            result_dir,
            run_manifest,
            errors,
        )
        if (
            structure_checks
            and confidence_checks
            and len(structure_checks) != len(confidence_checks)
        ):
            errors.append(
                "Protenix 结构文件数量与 summary confidence 文件数量不一致"
            )
        output_pairing = validate_output_pairing(
            structure_checks,
            confidence_checks,
            run_manifest,
            errors,
        )
    if status == "completed" and run_manifest.get("exit_code") != 0:
        errors.append("completed 要求 Protenix exit_code=0")
    if status == "completed" and run_manifest.get("failure") is not None:
        errors.append("completed 不得同时包含 failure")

    scanned_files = scan_result_safety(result_dir, errors)
    output_contract_ready = bool(structure_checks and confidence_checks) and not errors
    completion_ready = (
        run_manifest.get("inference_allowed") is True
        and run_manifest.get("exit_code") == 0
        and output_contract_ready
    )
    if status == "completed" and not completion_ready:
        errors.append("completed 未满足必需输出和 Validator 条件")

    if status == "completed" and not errors:
        review_disposition = (
            "review_required"
            if (
                ligand_handoff.get("review_disposition") == "review_required"
                or disposition == "review_required"
                or bool(scientific_review_reasons)
            )
            else "ready_for_review"
        )
    elif status == "not_run":
        review_disposition = (
            "rejected"
            if ligand_handoff.get("review_disposition") == "rejected"
            else "review_required"
        )
    else:
        review_disposition = "rejected" if errors else "review_required"

    file_checks = []
    for name in sorted(set(REQUIRED_JSON) | {"protenix_input.json"}):
        path = result_dir / name
        file_checks.append(
            {
                "file": name,
                "status": "present" if path.is_file() else "missing",
                "sha256": sha256_file(path) if path.is_file() else None,
            }
        )
    for item in (
        structure_checks
        + auxiliary_structure_checks
        + confidence_checks
        + full_confidence_checks
    ):
        file_checks.append(
            {
                "file": item.get("file"),
                "status": item.get("status"),
                "sha256": item.get("sha256"),
            }
        )

    report = {
        "schema_version": SCHEMA_VERSION,
        "workflow": WORKFLOW,
        "generated_at_utc": generated_at_utc or now_utc(),
        "job_id": input_manifest.get("job_id"),
        "status": status,
        "review_disposition": review_disposition,
        "valid": not errors,
        "completion_ready": completion_ready,
        "input_check": input_check,
        "structure_checks": structure_checks,
        "auxiliary_structure_checks": auxiliary_structure_checks,
        "confidence_checks": confidence_checks,
        "full_confidence_checks": full_confidence_checks,
        "output_pairing": output_pairing,
        "output_inventory": output_inventory,
        "file_checks": file_checks,
        "scanned_files": scanned_files,
        "errors": errors,
        "warnings": sorted(set(warnings)),
        "notices": notices,
        "scientific_boundaries": [
            "预测结构不是实验结构。",
            "ranking score、pLDDT、pTM、ipTM 不是 Kd、IC50、结合自由能、药效或实验准确率。",
            "高置信度不能排除错误构象、错误质子化、错误配体形式或错误生物学组装。",
            "专业数据集命中和摘要不等于论文全文已经核验。",
            "结果必须由结构生物学、计算化学或药物研发人员复核。",
            "本 Skill 不输出实验执行、安全、临床或治疗结论。",
        ],
    }
    report["result_fingerprint"] = document_fingerprint(report)
    if write_artifacts:
        write_json(result_dir / "validation_report.json", report)
        html_path = result_dir / "report.html"
        write_text_safely(html_path, render_report(report))
    return report


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="独立校验 Protenix 结果目录并生成专家复核报告。"
    )
    parser.add_argument("result_dir", type=Path)
    parser.add_argument("--generated-at")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        report = validate_result_dir(
            args.result_dir,
            generated_at_utc=args.generated_at,
            write_artifacts=True,
        )
    except (OSError, ValidationDependencyFailure) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 3
    print(json.dumps(report, ensure_ascii=False, indent=2))
    return 0 if report["valid"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
