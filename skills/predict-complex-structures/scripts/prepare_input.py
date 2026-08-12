#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import math
import re
import shutil
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


SCHEMA_VERSION = "1.0.0"
WORKFLOW = "predict-complex-structures"
PROTEIN_ALPHABET = frozenset("ACDEFGHIKLMNPQRSTVWYX")
DNA_ALPHABET = frozenset("ATGCN")
RNA_ALPHABET = frozenset("AUGCN")
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
JOB_ID_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]{0,127}$")
CHAIN_ID_RE = re.compile(r"^[A-Za-z0-9]{1,4}$")
CCD_RE = re.compile(r"^CCD_[A-Za-z0-9_]+$")
ION_RE = re.compile(r"^[A-Z][A-Z0-9]{0,7}$")
ATOM_RE = re.compile(r"^[A-Za-z0-9'_*+.-]{1,32}$")
SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
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
FILE_LIMIT_BYTES = 512 * 1024 * 1024
JSON_LIMIT_BYTES = 64 * 1024 * 1024
TEXT_SCIENTIFIC_SUFFIXES = {
    ".a3m",
    ".cif",
    ".fa",
    ".fasta",
    ".hhr",
    ".mmcif",
    ".mol",
    ".mol2",
    ".pdb",
    ".sdf",
}
SUPPORTED_MODELS = {
    "protenix-v2",
    "protenix_base_20250630_v1.0.0",
    "protenix_base_constraint_v0.5.0",
    "protenix_base_default_v0.5.0",
    "protenix_base_default_v1.0.0",
    "protenix_mini_default_v0.5.0",
    "protenix_tiny_default_v0.5.0",
}
TEMPLATE_AND_RNA_MSA_MODELS = {
    "protenix-v2",
    "protenix_base_20250630_v1.0.0",
    "protenix_base_default_v1.0.0",
}


class InputFailure(ValueError):
    pass


class DependencyFailure(RuntimeError):
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


def sha256_json(value: Any) -> str:
    return sha256_bytes(canonical_json(value).encode("utf-8"))


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
    return sha256_json(_without_nondeterministic_fields(document))


def serialized_safety_errors(
    value: Any,
    *,
    check_sensitive_paths: bool = True,
) -> list[str]:
    text = (
        value
        if isinstance(value, str)
        else json.dumps(value, ensure_ascii=False, sort_keys=True)
    )
    errors: list[str] = []
    if SECRET_RE.search(text):
        errors.append("possible credential detected")
    if check_sensitive_paths and SENSITIVE_PATH_RE.search(text):
        errors.append("host-sensitive path detected")
    return errors


def write_json(path: Path, document: Any) -> None:
    errors = serialized_safety_errors(document)
    if errors:
        raise InputFailure(f"{path.name}: {'; '.join(errors)}")
    temporary = path.with_name(f".{path.name}.tmp")
    temporary.write_text(
        json.dumps(document, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    temporary.replace(path)


def read_json(path: Path, label: str) -> tuple[Any, str]:
    if path.is_symlink() or not path.is_file():
        raise InputFailure(f"{label} 不是普通文件。")
    if path.stat().st_size > JSON_LIMIT_BYTES:
        raise InputFailure(f"{label} 超过 {JSON_LIMIT_BYTES} bytes。")
    raw = path.read_bytes()
    try:
        value = json.loads(raw.decode("utf-8"))
    except (UnicodeError, json.JSONDecodeError) as exc:
        raise InputFailure(f"{label} 不是有效 UTF-8 JSON：{exc}") from exc
    errors = serialized_safety_errors(value, check_sensitive_paths=False)
    if errors:
        raise InputFailure(f"{label}: {'; '.join(errors)}")
    return value, sha256_bytes(raw)


def resolve_input_path(base: Path, raw: Any, label: str) -> Path:
    if not isinstance(raw, str) or not raw.strip():
        raise InputFailure(f"{label} 必须是非空文件路径。")
    root = base.resolve()
    path = Path(raw)
    if not path.is_absolute():
        path = root / path
    if path.is_symlink():
        raise InputFailure(f"{label} 不允许 symlink。")
    resolved = path.resolve()
    try:
        resolved.relative_to(root)
    except ValueError as exc:
        raise InputFailure(f"{label} 必须位于 request 文件目录内。") from exc
    return resolved


def verify_expected_hash(actual: str, expected: Any, label: str) -> None:
    expected_text = str(expected).strip().lower()
    if not SHA256_RE.fullmatch(expected_text):
        raise InputFailure(f"{label}.sha256 必须提供 64 位小写十六进制。")
    if actual != expected_text:
        raise InputFailure(f"{label} SHA-256 与声明不一致。")


def validate_staged_file_safety(source: Path, label: str) -> None:
    if source.suffix.lower() not in TEXT_SCIENTIFIC_SUFFIXES:
        return
    try:
        with source.open(encoding="utf-8", errors="replace") as handle:
            for line_number, line in enumerate(handle, start=1):
                errors = serialized_safety_errors(line)
                if errors:
                    raise InputFailure(
                        f"{label}:{line_number}: {'; '.join(errors)}"
                    )
    except OSError as exc:
        raise InputFailure(f"{label} 无法执行安全扫描：{exc}") from exc


def stage_file(
    source: Path,
    result_dir: Path,
    category: str,
    expected_sha256: Any = None,
    allowed_suffixes: frozenset[str] | None = None,
) -> dict[str, Any]:
    if source.is_symlink() or not source.is_file():
        raise InputFailure(f"{category} 不是普通文件。")
    suffix = source.suffix.lower()
    if allowed_suffixes is not None and suffix not in allowed_suffixes:
        raise InputFailure(
            f"{category} 文件扩展名必须属于 {sorted(allowed_suffixes)}。"
        )
    size = source.stat().st_size
    if size <= 0 or size > FILE_LIMIT_BYTES:
        raise InputFailure(f"{category} 文件大小不在允许范围。")
    validate_staged_file_safety(source, category)
    actual_hash = sha256_file(source)
    verify_expected_hash(actual_hash, expected_sha256, category)
    safe_name = re.sub(r"[^A-Za-z0-9._-]+", "_", source.name).strip("._")
    if not safe_name:
        safe_name = f"{category}.bin"
    relative = Path("inputs") / f"{category}-{actual_hash[:12]}-{safe_name}"
    target = result_dir / relative
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(source, target)
    return {
        "role": category,
        "file_name": safe_name,
        "relative_path": relative.as_posix(),
        "sha256": actual_hash,
        "size_bytes": size,
    }


def normalize_sequence(value: Any, alphabet: frozenset[str], label: str) -> str:
    if not isinstance(value, str):
        raise InputFailure(f"{label} 必须是字符串。")
    sequence = "".join(value.split()).upper()
    if not sequence:
        raise InputFailure(f"{label} 不能为空。")
    invalid = sorted(set(sequence) - alphabet)
    if invalid:
        raise InputFailure(f"{label} 含非法字符：{''.join(invalid)}。")
    return sequence


def positive_integer(value: Any, label: str, maximum: int | None = None) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise InputFailure(f"{label} 必须是正整数。")
    if maximum is not None and value > maximum:
        raise InputFailure(f"{label} 超出实体范围。")
    return value


def finite_number(value: Any, label: str, *, minimum: float = 0.0) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise InputFailure(f"{label} 必须是数值。")
    number = float(value)
    if not math.isfinite(number) or number < minimum:
        raise InputFailure(f"{label} 必须是不小于 {minimum} 的有限数值。")
    return number


def normalize_modifications(
    raw: Any,
    *,
    entity_type: str,
    sequence_length: int,
    label: str,
) -> list[dict[str, Any]]:
    if raw is None:
        return []
    if not isinstance(raw, list):
        raise InputFailure(f"{label} 必须是 array。")
    if entity_type == "proteinChain":
        type_key, position_key = "ptmType", "ptmPosition"
    else:
        type_key, position_key = "modificationType", "basePosition"
    normalized: list[dict[str, Any]] = []
    for index, item in enumerate(raw):
        item_label = f"{label}[{index}]"
        if not isinstance(item, dict) or set(item) != {type_key, position_key}:
            raise InputFailure(
                f"{item_label} 只允许 {type_key}/{position_key}。"
            )
        code = str(item.get(type_key) or "").strip().upper()
        if not code.startswith("CCD_"):
            code = f"CCD_{code}"
        if not CCD_RE.fullmatch(code):
            raise InputFailure(f"{item_label}.{type_key} 不合法。")
        position = positive_integer(
            item.get(position_key),
            f"{item_label}.{position_key}",
            sequence_length,
        )
        normalized.append({type_key: code, position_key: position})
    return normalized


def expected_residue_names(
    sequence: str,
    residue_map: dict[str, str],
    modifications: list[dict[str, Any]],
    *,
    type_key: str,
    position_key: str,
) -> list[str]:
    """Build the residue identity contract, including declared CCD modifications."""
    names = [residue_map[residue] for residue in sequence]
    for modification in modifications:
        names[int(modification[position_key]) - 1] = str(
            modification[type_key]
        ).removeprefix("CCD_")
    # Protenix deliberately canonicalizes selenomethionine to methionine
    # before writing predictions (json_to_feature.py:mse_to_met).
    names = ["MET" if name == "MSE" else name for name in names]
    return names


def normalize_chain_id(value: Any, used: set[str], label: str) -> str:
    chain_id = str(value or "").strip()
    if not CHAIN_ID_RE.fullmatch(chain_id):
        raise InputFailure(f"{label} 必须是 1–4 位英文字母或数字。")
    if chain_id in used:
        raise InputFailure(f"链 ID 重复：{chain_id}。")
    used.add(chain_id)
    return chain_id


def load_rdkit() -> dict[str, Any]:
    try:
        import rdkit
        from rdkit import Chem
    except ImportError as exc:
        raise DependencyFailure(
            "缺少 rdkit==2025.9.2，无法建立配体校验口径。"
        ) from exc
    package_version = importlib.metadata.version("rdkit")
    if package_version != "2025.9.2":
        raise DependencyFailure(
            f"需要 rdkit==2025.9.2，当前为 {package_version}。"
        )
    return {"rdkit": rdkit, "Chem": Chem, "version": package_version}


def ligand_profile_from_molecule(molecule: Any) -> dict[str, Any]:
    heavy_elements = [
        atom.GetSymbol() for atom in molecule.GetAtoms() if atom.GetAtomicNum() > 1
    ]
    return {
        "heavy_atom_count": len(heavy_elements),
        "heavy_element_counts": dict(sorted(Counter(heavy_elements).items())),
        "formal_charge": sum(atom.GetFormalCharge() for atom in molecule.GetAtoms()),
    }


def ligand_profile_from_smiles(smiles: str, toolkit: dict[str, Any]) -> dict[str, Any]:
    molecule = toolkit["Chem"].MolFromSmiles(smiles)
    if molecule is None:
        raise InputFailure("标准化 Artifact 的配体结构无法被固定 RDKit 解析。")
    profile = ligand_profile_from_molecule(molecule)
    profile["canonical_isomeric_smiles"] = toolkit["Chem"].MolToSmiles(
        molecule,
        canonical=True,
        isomericSmiles=True,
    )
    profile["chiral_centers"] = [
        {"atom_index": index, "assignment": assignment}
        for index, assignment in toolkit["Chem"].FindMolChiralCenters(
            molecule,
            includeUnassigned=True,
            includeCIP=True,
        )
    ]
    return profile


def structure_identity(
    structure: Any,
    toolkit: dict[str, Any],
) -> dict[str, str] | None:
    if not isinstance(structure, str) or not structure.strip():
        return None
    molecule = toolkit["Chem"].MolFromSmiles(structure)
    if molecule is None and "M  END" in structure:
        molecule = toolkit["Chem"].MolFromMolBlock(structure, sanitize=True)
    if molecule is None:
        return None
    try:
        inchikey = toolkit["Chem"].MolToInchiKey(molecule)
    except (AttributeError, RuntimeError, ValueError):
        inchikey = ""
    return {
        "canonical_isomeric_smiles": toolkit["Chem"].MolToSmiles(
            molecule,
            canonical=True,
            isomericSmiles=True,
        ),
        "inchikey": inchikey,
    }


def ligand_profile_from_file(path: Path, toolkit: dict[str, Any]) -> dict[str, Any]:
    suffix = path.suffix.lower()
    Chem = toolkit["Chem"]
    if suffix in {".sdf", ".mol"}:
        molecule = Chem.MolFromMolFile(str(path), removeHs=False, sanitize=True)
    elif suffix == ".mol2":
        molecule = Chem.MolFromMol2File(str(path), removeHs=False, sanitize=True)
    elif suffix == ".pdb":
        molecule = Chem.MolFromPDBFile(str(path), removeHs=False, sanitize=True)
    else:
        raise InputFailure("结构文件只支持 SDF、MOL、MOL2 或 PDB。")
    if molecule is None:
        raise InputFailure("配体结构文件无法被固定 RDKit 解析。")
    if molecule.GetNumConformers() < 1 or not molecule.GetConformer().Is3D():
        raise InputFailure("FILE_ 配体结构必须包含 3D 构象。")
    profile = ligand_profile_from_molecule(molecule)
    profile["canonical_isomeric_smiles"] = Chem.MolToSmiles(
        molecule,
        canonical=True,
        isomericSmiles=True,
    )
    profile["chiral_centers"] = [
        {"atom_index": index, "assignment": assignment}
        for index, assignment in Chem.FindMolChiralCenters(
            molecule,
            includeUnassigned=True,
            includeCIP=True,
        )
    ]
    return profile


def verify_upstream_fingerprint(document: dict[str, Any], label: str) -> None:
    expected = document.get("result_fingerprint")
    if expected is None:
        raise InputFailure(f"{label} 缺少 result_fingerprint。")
    if expected != document_fingerprint(document):
        raise InputFailure(f"{label} result_fingerprint 与内容不一致。")


def find_identity_resolution(
    document: dict[str, Any],
    resolution_id: Any,
) -> dict[str, Any]:
    if document.get("workflow") != "chemical-identity-resolution":
        raise InputFailure("identity_artifact workflow 不正确。")
    verify_upstream_fingerprint(document, "identity_artifact")
    resolutions = document.get("resolutions")
    if not isinstance(resolutions, list) or not resolutions:
        raise InputFailure("identity_artifact.resolutions 不能为空。")
    if resolution_id is None and len(resolutions) == 1:
        return resolutions[0]
    for resolution in resolutions:
        request = resolution.get("request") if isinstance(resolution, dict) else None
        candidate_id = (
            request.get("id") if isinstance(request, dict) else resolution.get("id")
        )
        if candidate_id == resolution_id:
            return resolution
    raise InputFailure("identity_resolution_id 未在 identity_artifact 中找到。")


def find_standardization_record(
    document: dict[str, Any],
    record_id: Any,
) -> dict[str, Any]:
    if document.get("workflow") != "chemical-structure-standardization-qc":
        raise InputFailure("standardization_artifact workflow 不正确。")
    verify_upstream_fingerprint(document, "standardization_artifact")
    records = document.get("records")
    if not isinstance(records, list) or not records:
        raise InputFailure("standardization_artifact.records 不能为空。")
    if record_id is None and len(records) == 1:
        return records[0]
    for record in records:
        if isinstance(record, dict) and record.get("id") == record_id:
            return record
    raise InputFailure(
        "standardization_record_id 未在 standardization_artifact 中找到。"
    )


def load_evidence(
    request: dict[str, Any],
    request_dir: Path,
) -> tuple[dict[str, Any], dict[str, Any] | None]:
    raw_path = request.get("evidence_artifact")
    if raw_path is None:
        document = {
            "schema_version": SCHEMA_VERSION,
            "workflow": WORKFLOW,
            "evidence_status": "not_run",
            "source": "dataPro_search",
            "notice": "未提供专业数据集证据 Artifact；未执行或不可用。",
            "items": [],
        }
        document["evidence_fingerprint"] = document_fingerprint(document)
        return document, None

    path = resolve_input_path(request_dir, raw_path, "evidence_artifact")
    raw, file_hash = read_json(path, "evidence_artifact")
    if not isinstance(raw, dict):
        raise InputFailure("evidence_artifact 顶层必须是 object。")
    code = raw.get("code")
    status = (
        "completed"
        if code == 0 and raw.get("dataset_type") == "academic_search"
        else "error"
    )
    items = raw.get("items")
    if not isinstance(items, list):
        raise InputFailure("evidence_artifact.items 必须是 array。")
    document = {
        "schema_version": SCHEMA_VERSION,
        "workflow": WORKFLOW,
        "evidence_status": status,
        "source": "dataPro_search",
        "source_artifact_sha256": file_hash,
        "response": raw,
        "item_review": [
            {
                "index": index,
                "verification_status": "candidate_only",
                "doi_present": bool(
                    isinstance(item, dict)
                    and str(
                        item.get("doi")
                        or (
                            item.get("extra_data", {}).get("doi")
                            if isinstance(item.get("extra_data"), dict)
                            else ""
                        )
                    ).strip()
                ),
            }
            for index, item in enumerate(items)
        ],
        "notice": (
            "题名、摘要和搜索命中是候选证据；本 Skill 未核验论文全文。"
        ),
    }
    document["evidence_fingerprint"] = document_fingerprint(document)
    return document, {
        "role": "evidence_artifact",
        "file_name": path.name,
        "sha256": file_hash,
        "size_bytes": path.stat().st_size,
    }


def prepare_ligand(
    ligand: Any,
    request_dir: Path,
    result_dir: Path,
    used_chain_ids: set[str],
) -> tuple[
    dict[str, Any],
    dict[str, Any] | None,
    list[dict[str, Any]],
    bool,
    str,
]:
    if ligand is None:
        document = {
            "schema_version": SCHEMA_VERSION,
            "workflow": WORKFLOW,
            "ligand_status": "not_applicable",
            "review_disposition": "ready_for_review",
            "input_mode": None,
            "upstream": {},
            "selected_form": None,
            "expected_output_profile": None,
            "review_reasons": [],
        }
        document["ligand_handoff_fingerprint"] = document_fingerprint(document)
        return document, None, [], True, "ready_for_review"
    if not isinstance(ligand, dict):
        raise InputFailure("ligand 必须是 object。")

    chain_id = normalize_chain_id(
        ligand.get("chain_id", "L"),
        used_chain_ids,
        "ligand.chain_id",
    )
    input_mode = ligand.get("input_mode")
    if input_mode not in {"standardized_smiles", "ccd", "structure_file"}:
        raise InputFailure("ligand.input_mode 不受支持。")

    source_files: list[dict[str, Any]] = []
    upstream: dict[str, Any] = {}
    review_reasons: list[str] = []
    inference_allowed = True
    disposition = "ready_for_review"
    identity_handoff_record: dict[str, Any] | None = None

    identity_path_raw = ligand.get("identity_artifact")
    requires_identity = bool(
        str(ligand.get("reported_name") or "").strip()
        or (
            isinstance(ligand.get("source"), dict)
            and ligand["source"].get("identifier")
        )
    )
    if requires_identity and identity_path_raw is None and input_mode != "ccd":
        inference_allowed = False
        disposition = "review_required"
        review_reasons.append("R-LIGAND-IDENTITY-ARTIFACT-MISSING")
    if identity_path_raw is not None:
        identity_path = resolve_input_path(
            request_dir,
            identity_path_raw,
            "identity_artifact",
        )
        identity, identity_hash = read_json(identity_path, "identity_artifact")
        if not isinstance(identity, dict):
            raise InputFailure("identity_artifact 顶层必须是 object。")
        resolution = find_identity_resolution(
            identity,
            ligand.get("identity_resolution_id"),
        )
        handoff = resolution.get("standardization_handoff")
        handoff_records = (
            handoff.get("records") if isinstance(handoff, dict) else None
        )
        if (
            isinstance(handoff, dict)
            and handoff.get("status") == "ready"
            and isinstance(handoff_records, list)
            and len(handoff_records) == 1
            and isinstance(handoff_records[0], dict)
        ):
            identity_handoff_record = handoff_records[0]
        upstream["identity"] = {
            "artifact_sha256": identity_hash,
            "result_fingerprint": identity.get("result_fingerprint"),
            "resolution_id": ligand.get("identity_resolution_id"),
            "disposition": resolution.get("disposition"),
            "record_alignment_status": resolution.get("record_alignment_status"),
            "handoff_status": (
                handoff.get("status") if isinstance(handoff, dict) else None
            ),
            "handoff_record": (
                {
                    "id": identity_handoff_record.get("id"),
                    "structure": identity_handoff_record.get("structure"),
                    "source_inchikey": identity_handoff_record.get(
                        "source_inchikey"
                    ),
                }
                if identity_handoff_record is not None
                else None
            ),
        }
        source_files.append(
            {
                "role": "identity_artifact",
                "file_name": identity_path.name,
                "sha256": identity_hash,
                "size_bytes": identity_path.stat().st_size,
            }
        )
        if identity_handoff_record is None:
            inference_allowed = False
            disposition = "review_required"
            review_reasons.append("R-LIGAND-IDENTITY-NOT-CONFIRMED")

    standardization_record: dict[str, Any] | None = None
    standardization_path_raw = ligand.get("standardization_artifact")
    if input_mode in {"standardized_smiles", "structure_file"}:
        if standardization_path_raw is None:
            raise InputFailure(
                f"{input_mode} 必须提供 standardization_artifact。"
            )
        standardization_path = resolve_input_path(
            request_dir,
            standardization_path_raw,
            "standardization_artifact",
        )
        standardization, standardization_hash = read_json(
            standardization_path,
            "standardization_artifact",
        )
        if not isinstance(standardization, dict):
            raise InputFailure("standardization_artifact 顶层必须是 object。")
        standardization_record = find_standardization_record(
            standardization,
            ligand.get("standardization_record_id"),
        )
        upstream["standardization"] = {
            "artifact_sha256": standardization_hash,
            "result_fingerprint": standardization.get("result_fingerprint"),
            "record_id": standardization_record.get("id"),
            "profile": (
                standardization.get("options", {}).get("profile")
                if isinstance(standardization.get("options"), dict)
                else None
            ),
            "parse_status": standardization_record.get("parse_status"),
            "standardization_status": standardization_record.get(
                "standardization_status"
            ),
            "disposition": standardization_record.get("disposition"),
        }
        source_files.append(
            {
                "role": "standardization_artifact",
                "file_name": standardization_path.name,
                "sha256": standardization_hash,
                "size_bytes": standardization_path.stat().st_size,
            }
        )
        blocked = (
            standardization_record.get("disposition") == "rejected"
            or standardization_record.get("parse_status") == "error"
            or standardization_record.get("standardization_status")
            in {"not_run", "error"}
        )
        if blocked:
            inference_allowed = False
            disposition = "rejected"
            review_reasons.append("R-UPSTREAM-STANDARDIZATION-BLOCKED")
        elif standardization_record.get("disposition") == "review_required":
            disposition = "review_required"
            review_reasons.extend(
                str(item)
                for item in standardization_record.get(
                    "human_review_required",
                    [],
                )
            )
            review_reasons.append("R-UPSTREAM-REVIEW-REQUIRED")

    toolkit = load_rdkit()
    artifact_binding: dict[str, Any] = {
        "status": "not_applicable",
        "checks": {},
    }
    if (
        identity_handoff_record is not None
        and standardization_record is not None
    ):
        identity_structure = structure_identity(
            identity_handoff_record.get("structure"),
            toolkit,
        )
        original_structure = structure_identity(
            standardization_record.get("original_structure"),
            toolkit,
        )
        checks = {
            "record_id_match": (
                identity_handoff_record.get("id")
                == standardization_record.get("id")
            ),
            "canonical_structure_match": (
                identity_structure is not None
                and original_structure is not None
                and identity_structure["canonical_isomeric_smiles"]
                == original_structure["canonical_isomeric_smiles"]
            ),
            "source_inchikey_match": (
                identity_structure is not None
                and original_structure is not None
                and bool(identity_handoff_record.get("source_inchikey"))
                and identity_handoff_record.get("source_inchikey")
                == identity_structure["inchikey"]
                == original_structure["inchikey"]
            ),
        }
        artifact_binding = {
            "status": "matched" if all(checks.values()) else "mismatch",
            "checks": checks,
            "identity_record_id": identity_handoff_record.get("id"),
            "standardization_record_id": standardization_record.get("id"),
            "identity_structure": identity_structure,
            "standardization_original_structure": original_structure,
        }
        if artifact_binding["status"] != "matched":
            inference_allowed = False
            disposition = "rejected"
            review_reasons.append("R-UPSTREAM-ARTIFACT-MISMATCH")
    official_ligand: dict[str, Any]
    selected_form: dict[str, Any]
    output_profile: dict[str, Any] | None

    if input_mode == "standardized_smiles":
        assert standardization_record is not None
        smiles = standardization_record.get("standardized_structure")
        if not isinstance(smiles, str) or not smiles:
            inference_allowed = False
            disposition = "rejected"
            review_reasons.append("R-STANDARDIZED-STRUCTURE-MISSING")
            smiles = ""
            output_profile = None
        else:
            output_profile = ligand_profile_from_smiles(smiles, toolkit)
        selected_form = {
            "reported_form": standardization_record.get("original_structure"),
            "standardized_form": smiles or None,
            "parent_form": standardization_record.get("parent_structure"),
            "inchikey": standardization_record.get("inchikey"),
            "parent_inchikey": standardization_record.get("parent_inchikey"),
            "fragment_analysis": standardization_record.get(
                "fragment_analysis"
            ),
            "transformations": standardization_record.get("transformations"),
            "qc_findings": standardization_record.get("qc_findings"),
            "upstream_review_reasons": standardization_record.get(
                "human_review_required",
                [],
            ),
            "selection": "standardized_form",
            "parent_used_for_prediction": False,
        }
        official_ligand = {
            "ligand": smiles,
            "count": 1,
            "id": [chain_id],
        }
    elif input_mode == "ccd":
        ccd_code = str(ligand.get("ccd_id") or "").strip().upper()
        if not ccd_code.startswith("CCD_"):
            ccd_code = f"CCD_{ccd_code}"
        if not CCD_RE.fullmatch(ccd_code):
            raise InputFailure("ligand.ccd_id 不符合 CCD_<code>。")
        if ligand.get("confirmation_status") not in {
            "user_confirmed",
            "expert_confirmed",
        }:
            inference_allowed = False
            disposition = "review_required"
            review_reasons.append("R-CCD-NOT-CONFIRMED")
        selected_form = {
            "reported_form": ligand.get("reported_name"),
            "standardized_form": None,
            "parent_form": None,
            "ccd_id": ccd_code,
            "selection": "confirmed_ccd",
            "parent_used_for_prediction": False,
        }
        output_profile = {"ccd_id": ccd_code}
        official_ligand = {
            "ligand": ccd_code,
            "count": 1,
            "id": [chain_id],
        }
    else:
        assert standardization_record is not None
        source_path = resolve_input_path(
            request_dir,
            ligand.get("structure_file"),
            "ligand.structure_file",
        )
        staged = stage_file(
            source_path,
            result_dir,
            "ligand",
            ligand.get("sha256"),
            frozenset({".mol", ".mol2", ".pdb", ".sdf"}),
        )
        source_files.append(staged)
        staged_path = result_dir / staged["relative_path"]
        output_profile = ligand_profile_from_file(staged_path, toolkit)
        selected_form = {
            "reported_form": standardization_record.get("original_structure"),
            "standardized_form": standardization_record.get(
                "standardized_structure"
            ),
            "parent_form": standardization_record.get("parent_structure"),
            "inchikey": standardization_record.get("inchikey"),
            "parent_inchikey": standardization_record.get("parent_inchikey"),
            "fragment_analysis": standardization_record.get(
                "fragment_analysis"
            ),
            "transformations": standardization_record.get("transformations"),
            "qc_findings": standardization_record.get("qc_findings"),
            "upstream_review_reasons": standardization_record.get(
                "human_review_required",
                [],
            ),
            "selection": "confirmed_3d_structure_file",
            "structure_file_sha256": staged["sha256"],
            "parent_used_for_prediction": False,
        }
        official_ligand = {
            "ligand": f"FILE_{staged['relative_path']}",
            "count": 1,
            "id": [chain_id],
        }

    document = {
        "schema_version": SCHEMA_VERSION,
        "workflow": WORKFLOW,
        "ligand_status": (
            "blocked"
            if not inference_allowed
            else "review_required"
            if disposition == "review_required"
            else "ready"
        ),
        "review_disposition": disposition,
        "chain_id": chain_id,
        "input_mode": input_mode,
        "reported_name": ligand.get("reported_name"),
        "source": ligand.get("source"),
        "upstream": upstream,
        "artifact_binding": artifact_binding,
        "selected_form": selected_form,
        "expected_output_profile": output_profile,
        "review_reasons": sorted(set(review_reasons)),
        "notice": (
            "reported、standardized 与 parent 均保留；模型只使用 selection "
            "明确指定的具体形式。parent 未用于预测。"
        ),
    }
    document["ligand_handoff_fingerprint"] = document_fingerprint(document)
    official_entity = {"ligand": official_ligand} if inference_allowed else None
    return document, official_entity, source_files, inference_allowed, disposition


def stage_chain_feature(
    chain: dict[str, Any],
    field: str,
    official_field: str,
    request_dir: Path,
    result_dir: Path,
    source_files: list[dict[str, Any]],
    role: str,
) -> tuple[str, str] | None:
    value = chain.get(field)
    if value is None:
        return None
    if not isinstance(value, dict):
        raise InputFailure(f"{field} 必须包含 path 和 sha256。")
    path = resolve_input_path(request_dir, value.get("path"), field)
    allowed_suffixes = (
        frozenset({".a3m", ".hhr"})
        if field == "template"
        else frozenset({".a3m"})
    )
    staged = stage_file(
        path,
        result_dir,
        role,
        value.get("sha256"),
        allowed_suffixes,
    )
    source_files.append(staged)
    return official_field, staged["relative_path"]


def build_entities(
    request: dict[str, Any],
    request_dir: Path,
    result_dir: Path,
    used_chain_ids: set[str],
    source_files: list[dict[str, Any]],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    proteins = request.get("proteins")
    if not isinstance(proteins, list) or not proteins:
        raise InputFailure("proteins 必须是非空 array。")
    entities: list[dict[str, Any]] = []
    polymer_manifest: list[dict[str, Any]] = []
    for index, item in enumerate(proteins):
        if not isinstance(item, dict):
            raise InputFailure(f"proteins[{index}] 必须是 object。")
        chain_id = normalize_chain_id(
            item.get("chain_id"),
            used_chain_ids,
            f"proteins[{index}].chain_id",
        )
        sequence = normalize_sequence(
            item.get("sequence"),
            PROTEIN_ALPHABET,
            f"proteins[{index}].sequence",
        )
        source = item.get("source")
        if not isinstance(source, dict) or not source.get("database"):
            raise InputFailure(f"proteins[{index}].source.database 必填。")
        organism = str(item.get("organism") or "").strip()
        if not organism:
            raise InputFailure(f"proteins[{index}].organism 必填。")
        chain = {
            "sequence": sequence,
            "count": 1,
            "id": [chain_id],
            "modifications": normalize_modifications(
                item.get("modifications"),
                entity_type="proteinChain",
                sequence_length=len(sequence),
                label=f"proteins[{index}].modifications",
            ),
        }
        for staged in (
            stage_chain_feature(
                item,
                "paired_msa",
                "pairedMsaPath",
                request_dir,
                result_dir,
                source_files,
                f"protein-{chain_id}-paired-msa",
            ),
            stage_chain_feature(
                item,
                "unpaired_msa",
                "unpairedMsaPath",
                request_dir,
                result_dir,
                source_files,
                f"protein-{chain_id}-unpaired-msa",
            ),
            stage_chain_feature(
                item,
                "template",
                "templatesPath",
                request_dir,
                result_dir,
                source_files,
                f"protein-{chain_id}-template",
            ),
        ):
            if staged is not None:
                chain[staged[0]] = staged[1]
        entities.append({"proteinChain": chain})
        polymer_manifest.append(
            {
                "chain_id": chain_id,
                "entity_type": "proteinChain",
                "sequence": sequence,
                "sequence_sha256": sha256_bytes(sequence.encode("ascii")),
                "length": len(sequence),
                "expected_residue_names": expected_residue_names(
                    sequence,
                    PROTEIN_RESIDUE_NAMES,
                    chain["modifications"],
                    type_key="ptmType",
                    position_key="ptmPosition",
                ),
                "organism": organism,
                "source": source,
                "modifications": chain["modifications"],
            }
        )

    for field, official_key, alphabet in (
        ("dna", "dnaSequence", DNA_ALPHABET),
        ("rna", "rnaSequence", RNA_ALPHABET),
    ):
        raw_entities = request.get(field) or []
        if not isinstance(raw_entities, list):
            raise InputFailure(f"{field} 必须是 array。")
        for index, item in enumerate(raw_entities):
            if not isinstance(item, dict):
                raise InputFailure(f"{field}[{index}] 必须是 object。")
            chain_id = normalize_chain_id(
                item.get("chain_id"),
                used_chain_ids,
                f"{field}[{index}].chain_id",
            )
            sequence = normalize_sequence(
                item.get("sequence"),
                alphabet,
                f"{field}[{index}].sequence",
            )
            official = {
                "sequence": sequence,
                "count": 1,
                "id": [chain_id],
                "modifications": normalize_modifications(
                    item.get("modifications"),
                    entity_type=official_key,
                    sequence_length=len(sequence),
                    label=f"{field}[{index}].modifications",
                ),
            }
            if field == "rna":
                staged = stage_chain_feature(
                    item,
                    "unpaired_msa",
                    "unpairedMsaPath",
                    request_dir,
                    result_dir,
                    source_files,
                    f"rna-{chain_id}-unpaired-msa",
                )
                if staged is not None:
                    official[staged[0]] = staged[1]
            entities.append({official_key: official})
            source = item.get("source")
            if not isinstance(source, dict) or not source.get("database"):
                raise InputFailure(f"{field}[{index}].source.database 必填。")
            organism = str(item.get("organism") or "").strip()
            if not organism:
                raise InputFailure(f"{field}[{index}].organism 必填。")
            polymer_manifest.append(
                {
                    "chain_id": chain_id,
                    "entity_type": official_key,
                    "sequence": sequence,
                    "sequence_sha256": sha256_bytes(sequence.encode("ascii")),
                    "length": len(sequence),
                    "organism": organism,
                    "source": source,
                    "modifications": official["modifications"],
                    "expected_residue_names": expected_residue_names(
                        sequence,
                        DNA_RESIDUE_NAMES if field == "dna" else RNA_RESIDUE_NAMES,
                        official["modifications"],
                        type_key="modificationType",
                        position_key="basePosition",
                    ),
                }
            )

    ions = request.get("ions") or []
    if not isinstance(ions, list):
        raise InputFailure("ions 必须是 array。")
    for index, item in enumerate(ions):
        if not isinstance(item, dict):
            raise InputFailure(f"ions[{index}] 必须是 object。")
        chain_id = normalize_chain_id(
            item.get("chain_id"),
            used_chain_ids,
            f"ions[{index}].chain_id",
        )
        ion = str(item.get("ion") or "").strip().upper()
        if not ION_RE.fullmatch(ion) or ion.startswith("CCD_"):
            raise InputFailure(f"ions[{index}].ion 不符合官方离子 code。")
        entities.append({"ion": {"ion": ion, "count": 1, "id": [chain_id]}})

    cofactors = request.get("cofactors") or []
    if not isinstance(cofactors, list):
        raise InputFailure("cofactors 必须是 array。")
    for index, item in enumerate(cofactors):
        if not isinstance(item, dict):
            raise InputFailure(f"cofactors[{index}] 必须是 object。")
        chain_id = normalize_chain_id(
            item.get("chain_id"),
            used_chain_ids,
            f"cofactors[{index}].chain_id",
        )
        ccd = str(item.get("ccd_id") or "").strip().upper()
        if not ccd.startswith("CCD_"):
            ccd = f"CCD_{ccd}"
        if not CCD_RE.fullmatch(ccd):
            raise InputFailure(f"cofactors[{index}].ccd_id 不合法。")
        entities.append(
            {"ligand": {"ligand": ccd, "count": 1, "id": [chain_id]}}
        )
    return entities, polymer_manifest


def normalize_model(raw: Any) -> dict[str, Any]:
    if not isinstance(raw, dict):
        raise InputFailure("model 必须是 object。")
    name = str(raw.get("name") or "").strip()
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]{0,127}", name):
        raise InputFailure("model.name 不合法。")
    if name not in SUPPORTED_MODELS:
        raise InputFailure(
            "model.name 不在首版受控模型清单；ESM/ISM 等额外权重模型尚未支持。"
        )
    seeds = raw.get("seeds")
    if seeds is None and raw.get("seed") is not None:
        seeds = [raw["seed"]]
    if not isinstance(seeds, list) or not seeds:
        raise InputFailure("model.seeds 必须是非空整数 array。")
    normalized_seeds: list[int] = []
    for seed in seeds:
        if (
            isinstance(seed, bool)
            or not isinstance(seed, int)
            or seed < 0
            or seed > 2**32 - 1
        ):
            raise InputFailure("model.seeds 只能包含 0–2^32-1 的整数。")
        normalized_seeds.append(seed)
    if len(set(normalized_seeds)) != len(normalized_seeds):
        raise InputFailure("model.seeds 不得重复。")
    for key in ("msa_strategy", "template_strategy", "rna_msa_strategy"):
        if raw.get(key, "none") not in {"none", "precomputed"}:
            raise InputFailure(f"model.{key} 只允许 none/precomputed。")
    if (
        raw.get("template_strategy", "none") != "none"
        and name not in TEMPLATE_AND_RNA_MSA_MODELS
    ):
        raise InputFailure(f"{name} 不支持模板推理。")
    if (
        raw.get("rna_msa_strategy", "none") != "none"
        and name not in TEMPLATE_AND_RNA_MSA_MODELS
    ):
        raise InputFailure(f"{name} 不支持 RNA MSA 推理。")
    dtype = raw.get("dtype", "bf16")
    if dtype not in {"bf16", "fp32"}:
        raise InputFailure("model.dtype 只允许 bf16/fp32。")
    numeric: dict[str, int] = {}
    for key, default, maximum in (
        ("samples", 5, 100),
        ("steps", 200, 10_000),
        ("cycles", 10, 100),
    ):
        value = raw.get(key, default)
        if (
            isinstance(value, bool)
            or not isinstance(value, int)
            or value <= 0
            or value > maximum
        ):
            raise InputFailure(f"model.{key} 必须是 1–{maximum} 的整数。")
        numeric[key] = value
    weights_id = str(raw.get("weights_id") or name).strip()
    if weights_id != name:
        raise InputFailure(
            "model.weights_id 必须与 model.name 一致；当前 Runner 只通过 "
            "Protenix --model_name 选择官方权重，不得记录未实际使用的权重标识。"
        )
    return {
        "name": name,
        "version": str(raw.get("version") or "not_exposed_by_cli"),
        "weights_id": weights_id,
        "seeds": normalized_seeds,
        "msa_strategy": raw.get("msa_strategy", "none"),
        "template_strategy": raw.get("template_strategy", "none"),
        "rna_msa_strategy": raw.get("rna_msa_strategy", "none"),
        "dtype": dtype,
        **numeric,
    }


def validate_feature_strategies(
    model: dict[str, Any],
    entities: list[dict[str, Any]],
) -> None:
    proteins = [
        entity["proteinChain"]
        for entity in entities
        if "proteinChain" in entity
    ]
    rna_chains = [
        entity["rnaSequence"] for entity in entities if "rnaSequence" in entity
    ]
    protein_msa_present = [
        "pairedMsaPath" in chain or "unpairedMsaPath" in chain
        for chain in proteins
    ]
    templates_present = ["templatesPath" in chain for chain in proteins]
    rna_msa_present = ["unpairedMsaPath" in chain for chain in rna_chains]

    if model["msa_strategy"] == "precomputed" and not all(
        protein_msa_present
    ):
        raise InputFailure(
            "msa_strategy=precomputed 要求每条蛋白链提供 paired_msa 或 "
            "unpaired_msa。"
        )
    if model["msa_strategy"] == "none" and any(protein_msa_present):
        raise InputFailure("提供了蛋白 MSA 文件，但 msa_strategy=none。")
    if model["template_strategy"] == "precomputed" and not all(
        templates_present
    ):
        raise InputFailure(
            "template_strategy=precomputed 要求每条蛋白链提供 template。"
        )
    if model["template_strategy"] == "none" and any(templates_present):
        raise InputFailure("提供了模板文件，但 template_strategy=none。")
    if model["rna_msa_strategy"] == "precomputed":
        if not rna_chains or not all(rna_msa_present):
            raise InputFailure(
                "rna_msa_strategy=precomputed 要求每条 RNA 链提供 "
                "unpaired_msa。"
            )
    if model["rna_msa_strategy"] == "none" and any(rna_msa_present):
        raise InputFailure("提供了 RNA MSA 文件，但 rna_msa_strategy=none。")


def validate_covalent_bonds(
    raw: Any,
    entities: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    if raw is None:
        return []
    if not isinstance(raw, list):
        raise InputFailure("covalent_bonds 必须是 array。")
    normalized: list[dict[str, Any]] = []
    allowed = {
        "entity1",
        "copy1",
        "position1",
        "atom1",
        "entity2",
        "copy2",
        "position2",
        "atom2",
    }
    for index, item in enumerate(raw):
        label = f"covalent_bonds[{index}]"
        if not isinstance(item, dict) or set(item) - allowed:
            raise InputFailure(f"{label} 含非官方字段。")
        if not {
            "entity1",
            "position1",
            "atom1",
            "entity2",
            "position2",
            "atom2",
        } <= set(item):
            raise InputFailure(f"{label} 缺少定位原子的必需字段。")
        copy_keys = {"copy1", "copy2"} & set(item)
        if copy_keys not in (set(), {"copy1", "copy2"}):
            raise InputFailure(f"{label}.copy1/copy2 必须同时提供。")
        output: dict[str, Any] = {}
        for side in (1, 2):
            entity_number = positive_integer(
                item.get(f"entity{side}"),
                f"{label}.entity{side}",
                len(entities),
            )
            entity = entities[entity_number - 1]
            payload = next(iter(entity.values()))
            count = int(payload.get("count", 1))
            position_max = (
                len(payload.get("sequence", ""))
                if isinstance(payload.get("sequence"), str)
                else 1
            )
            output[f"entity{side}"] = entity_number
            if copy_keys:
                output[f"copy{side}"] = positive_integer(
                    item.get(f"copy{side}"),
                    f"{label}.copy{side}",
                    count,
                )
            output[f"position{side}"] = str(
                positive_integer(
                    int(item.get(f"position{side}"))
                    if isinstance(item.get(f"position{side}"), str)
                    and item.get(f"position{side}").isdigit()
                    else item.get(f"position{side}"),
                    f"{label}.position{side}",
                    position_max,
                )
            )
            atom = str(item.get(f"atom{side}") or "").strip()
            if not ATOM_RE.fullmatch(atom):
                raise InputFailure(f"{label}.atom{side} 不合法。")
            output[f"atom{side}"] = atom
        if output["entity1"] == output["entity2"] and output.get(
            "copy1", 1
        ) == output.get("copy2", 1) and output["position1"] == output[
            "position2"
        ] and output["atom1"] == output["atom2"]:
            raise InputFailure(f"{label} 不得将原子与自身成键。")
        normalized.append(output)
    return normalized


def validate_constraint(
    raw: Any,
    entities: list[dict[str, Any]],
) -> dict[str, Any] | None:
    if raw is None:
        return None
    if not isinstance(raw, dict) or not raw or set(raw) - {"contact", "pocket"}:
        raise InputFailure("constraint 只允许非空 contact/pocket。")

    def validate_locator(
        item: Any,
        label: str,
        *,
        atom_optional: bool,
    ) -> dict[str, Any]:
        if not isinstance(item, dict):
            raise InputFailure(f"{label} 必须是 object。")
        required = {"entity", "copy", "position"}
        allowed = required | ({"atom"} if atom_optional else set())
        if not required <= set(item) or set(item) - allowed:
            raise InputFailure(f"{label} 定位字段不合法。")
        entity_number = positive_integer(
            item.get("entity"),
            f"{label}.entity",
            len(entities),
        )
        payload = next(iter(entities[entity_number - 1].values()))
        count = int(payload.get("count", 1))
        position_max = (
            len(payload.get("sequence", ""))
            if isinstance(payload.get("sequence"), str)
            else 1
        )
        output = {
            "entity": entity_number,
            "copy": positive_integer(item.get("copy"), f"{label}.copy", count),
            "position": positive_integer(
                item.get("position"),
                f"{label}.position",
                position_max,
            ),
        }
        if "atom" in item:
            atom = str(item["atom"] or "").strip()
            if not ATOM_RE.fullmatch(atom):
                raise InputFailure(f"{label}.atom 不合法。")
            output["atom"] = atom
        return output

    normalized: dict[str, Any] = {}
    if "contact" in raw:
        contacts = raw["contact"]
        if not isinstance(contacts, list) or not contacts:
            raise InputFailure("constraint.contact 必须是非空 array。")
        normalized_contacts = []
        for index, item in enumerate(contacts):
            label = f"constraint.contact[{index}]"
            if not isinstance(item, dict):
                raise InputFailure(f"{label} 必须是 object。")
            side1 = {
                key.removesuffix("1"): value
                for key, value in item.items()
                if key in {"entity1", "copy1", "position1", "atom1"}
            }
            side2 = {
                key.removesuffix("2"): value
                for key, value in item.items()
                if key in {"entity2", "copy2", "position2", "atom2"}
            }
            left = validate_locator(side1, f"{label}.side1", atom_optional=True)
            right = validate_locator(side2, f"{label}.side2", atom_optional=True)
            max_distance = finite_number(
                item.get("max_distance"),
                f"{label}.max_distance",
            )
            min_distance = finite_number(
                item.get("min_distance", 0),
                f"{label}.min_distance",
            )
            if min_distance > max_distance:
                raise InputFailure(f"{label}: min_distance 不得大于 max_distance。")
            output = {
                "entity1": left["entity"],
                "copy1": left["copy"],
                "position1": left["position"],
                "entity2": right["entity"],
                "copy2": right["copy"],
                "position2": right["position"],
                "max_distance": max_distance,
                "min_distance": min_distance,
            }
            if "atom" in left:
                output["atom1"] = left["atom"]
            if "atom" in right:
                output["atom2"] = right["atom"]
            allowed_keys = set(output)
            if set(item) - allowed_keys:
                raise InputFailure(f"{label} 含非官方字段。")
            normalized_contacts.append(output)
        normalized["contact"] = normalized_contacts
    if "pocket" in raw:
        pocket = raw["pocket"]
        if not isinstance(pocket, dict) or set(pocket) != {
            "binder_chain",
            "contact_residues",
            "max_distance",
        }:
            raise InputFailure("constraint.pocket 字段不合法。")
        binder = pocket["binder_chain"]
        if not isinstance(binder, dict) or set(binder) != {"entity", "copy"}:
            raise InputFailure("constraint.pocket.binder_chain 字段不合法。")
        entity_number = positive_integer(
            binder.get("entity"),
            "constraint.pocket.binder_chain.entity",
            len(entities),
        )
        payload = next(iter(entities[entity_number - 1].values()))
        normalized_binder = {
            "entity": entity_number,
            "copy": positive_integer(
                binder.get("copy"),
                "constraint.pocket.binder_chain.copy",
                int(payload.get("count", 1)),
            ),
        }
        residues = pocket["contact_residues"]
        if not isinstance(residues, list) or not residues:
            raise InputFailure("constraint.pocket.contact_residues 必须是非空 array。")
        normalized["pocket"] = {
            "binder_chain": normalized_binder,
            "contact_residues": [
                validate_locator(
                    item,
                    f"constraint.pocket.contact_residues[{index}]",
                    atom_optional=False,
                )
                for index, item in enumerate(residues)
            ],
            "max_distance": finite_number(
                pocket["max_distance"],
                "constraint.pocket.max_distance",
            ),
        }
    return normalized


def prepare_job(
    request: dict[str, Any],
    request_dir: Path,
    result_dir: Path,
    generated_at_utc: str | None = None,
) -> tuple[dict[str, Any], int]:
    if not isinstance(request, dict):
        raise InputFailure("request 顶层必须是 object。")
    if request.get("schema_version") != SCHEMA_VERSION:
        raise InputFailure(f"schema_version 必须是 {SCHEMA_VERSION}。")
    job_id = str(request.get("job_id") or "").strip()
    if not JOB_ID_RE.fullmatch(job_id):
        raise InputFailure("job_id 不合法。")
    research_question = str(request.get("research_question") or "").strip()
    if not research_question:
        raise InputFailure("research_question 不能为空。")
    model = normalize_model(request.get("model"))

    if result_dir.is_symlink():
        raise InputFailure("result_dir 不允许 symlink。")
    if result_dir.exists() and any(result_dir.iterdir()):
        raise InputFailure("result_dir 必须不存在或为空。")
    result_dir.mkdir(parents=True, exist_ok=True)

    used_chain_ids: set[str] = set()
    source_files: list[dict[str, Any]] = []
    entities, polymers = build_entities(
        request,
        request_dir,
        result_dir,
        used_chain_ids,
        source_files,
    )
    ligand_handoff, ligand_entity, ligand_files, allowed, disposition = (
        prepare_ligand(
            request.get("ligand"),
            request_dir,
            result_dir,
            used_chain_ids,
        )
    )
    source_files.extend(ligand_files)
    if ligand_entity is not None:
        entities.append(ligand_entity)
    validate_feature_strategies(model, entities)

    reference_manifest: list[dict[str, Any]] = []
    references = request.get("reference_structures") or []
    if not isinstance(references, list):
        raise InputFailure("reference_structures 必须是 array。")
    for index, item in enumerate(references):
        if not isinstance(item, dict):
            raise InputFailure(f"reference_structures[{index}] 必须是 object。")
        source = resolve_input_path(
            request_dir,
            item.get("path"),
            f"reference_structures[{index}].path",
        )
        staged = stage_file(
            source,
            result_dir,
            f"reference-{index + 1}",
            item.get("sha256"),
            frozenset({".cif", ".mmcif", ".pdb"}),
        )
        source_files.append(staged)
        reference_manifest.append(
            {
                **staged,
                "source": item.get("source"),
                "identifier": item.get("identifier"),
                "use": "reference_or_validation_only",
            }
        )

    evidence, evidence_file = load_evidence(request, request_dir)
    if evidence_file is not None:
        source_files.append(evidence_file)

    covalent_bonds = validate_covalent_bonds(
        request.get("covalent_bonds"),
        entities,
    )
    constraint = validate_constraint(request.get("constraint"), entities)
    if constraint is not None:
        if "constraint" not in model["name"]:
            raise InputFailure(
                "当前 model 未声明 constraint 支持；不得静默传入软约束。"
            )

    scientific_review_reasons: list[str] = []
    if covalent_bonds:
        scientific_review_reasons.append(
            "COVALENT_BOND_NOT_INDEPENDENTLY_VERIFIED"
        )
    if constraint is not None:
        scientific_review_reasons.append(
            "SOFT_CONSTRAINT_NOT_INDEPENDENTLY_VERIFIED"
        )
    if allowed and scientific_review_reasons:
        disposition = "review_required"

    protenix_task: dict[str, Any] = {
        "name": job_id,
        "sequences": entities,
        "covalent_bonds": covalent_bonds,
    }
    if constraint is not None:
        protenix_task["constraint"] = constraint
    protenix_input = [protenix_task] if allowed else []
    protenix_input_sha256 = sha256_json(protenix_input)

    generated_at = generated_at_utc or now_utc()
    input_manifest = {
        "schema_version": SCHEMA_VERSION,
        "workflow": WORKFLOW,
        "generated_at_utc": generated_at,
        "job_id": job_id,
        "research_question": research_question,
        "polymers": polymers,
        "entity_chain_ids": sorted(used_chain_ids),
        "reference_structures": reference_manifest,
        "source_files": sorted(
            source_files,
            key=lambda item: (str(item.get("role")), str(item.get("sha256"))),
        ),
        "model_request": model,
        "covalent_bonds": covalent_bonds,
        "constraint": constraint,
        "scientific_review_reasons": scientific_review_reasons,
        "protenix_input_sha256": protenix_input_sha256,
        "bound_evidence_fingerprint": evidence["evidence_fingerprint"],
        "bound_ligand_handoff_fingerprint": ligand_handoff[
            "ligand_handoff_fingerprint"
        ],
    }
    input_manifest["input_fingerprint"] = document_fingerprint(input_manifest)

    run_manifest = {
        "schema_version": SCHEMA_VERSION,
        "workflow": WORKFLOW,
        "generated_at_utc": generated_at,
        "job_id": job_id,
        "status": "not_run",
        "review_disposition": disposition,
        "scientific_review_reasons": scientific_review_reasons,
        "inference_allowed": allowed,
        "input_fingerprint": input_manifest["input_fingerprint"],
        "model": model,
        "runtime": {
            "backend": "protenix_cli",
            "environment": "independent_gpu_runtime",
            "protenix_version": None,
            "device": None,
        },
        "command": None,
        "started_at_utc": None,
        "finished_at_utc": None,
        "elapsed_seconds": None,
        "exit_code": None,
        "failure": (
            None
            if allowed
            else {
                "code": "UPSTREAM_BLOCKED",
                "message": "配体身份或结构上游状态阻止 Protenix 推理。",
            }
        ),
        "outputs": {
            "structure_files": [],
            "auxiliary_structure_files": [],
            "summary_confidence_files": [],
            "full_confidence_files": [],
        },
        "notice": (
            "输入准备完成不表示 Protenix 已运行；只有独立 Validator 通过后才可 "
            "标记 completed。"
        ),
    }
    run_manifest["run_fingerprint"] = document_fingerprint(run_manifest)

    for path, document in (
        (result_dir / "input_manifest.json", input_manifest),
        (result_dir / "evidence.json", evidence),
        (result_dir / "ligand_handoff.json", ligand_handoff),
        (result_dir / "protenix_input.json", protenix_input),
        (result_dir / "run_manifest.json", run_manifest),
    ):
        write_json(path, document)

    result = {
        "input_manifest": input_manifest,
        "evidence": evidence,
        "ligand_handoff": ligand_handoff,
        "protenix_input": protenix_input,
        "run_manifest": run_manifest,
    }
    return result, 0 if allowed else 2


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="准备可审计的官方兼容 Protenix 输入和运行清单。"
    )
    parser.add_argument("--request", required=True, type=Path)
    parser.add_argument("--result-dir", required=True, type=Path)
    parser.add_argument("--generated-at")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        if args.request.is_symlink():
            raise InputFailure("request 不允许 symlink。")
        request_path = args.request.absolute()
        request, _ = read_json(request_path, "request")
        if not isinstance(request, dict):
            raise InputFailure("request 顶层必须是 object。")
        _, exit_code = prepare_job(
            request,
            request_path.parent,
            args.result_dir.absolute(),
            generated_at_utc=args.generated_at,
        )
    except (InputFailure, DependencyFailure, OSError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 3
    print(f"prepared Protenix result package: {args.result_dir}")
    return exit_code


if __name__ == "__main__":
    raise SystemExit(main())
