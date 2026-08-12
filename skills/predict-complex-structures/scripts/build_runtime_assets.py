#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
import sys
from pathlib import Path
from typing import Any


SCHEMA_VERSION = "1.0.0"
MODEL_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]{0,127}$")
BASE_CACHES = {
    "ccd_components_file": "common/components.cif",
    "ccd_components_rdkit_mol_file": "common/components.cif.rdkit_mol.pkl",
    "obsolete_release_data_csv": "common/obsolete_release_date.csv",
    "pdb_cluster_file": "common/clusters-by-entity-40.txt",
}
TEMPLATE_CACHES = {
    "obsolete_pdbs_path": "common/obsolete_to_successor.json",
    "release_dates_path": "common/release_date_cache.json",
}


class AssetFailure(RuntimeError):
    pass


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def checked_file(root: Path, relative: str, label: str) -> dict[str, Any]:
    unresolved_path = root / relative
    if unresolved_path.is_symlink():
        raise AssetFailure(f"{label} 不允许 symlink：{relative}。")
    path = unresolved_path.resolve()
    try:
        path.relative_to(root)
    except ValueError as exc:
        raise AssetFailure(f"{label} 路径逃逸 Protenix 根目录。") from exc
    if not path.is_file() or path.stat().st_size <= 0:
        raise AssetFailure(f"{label} 文件缺失、为空或为 symlink：{relative}。")
    return {
        "name": label,
        "relative_path": relative,
        "sha256": sha256_file(path),
    }


def detect_version(command: str) -> str:
    try:
        completed = subprocess.run(
            [command, "--version"],
            capture_output=True,
            text=True,
            check=False,
            timeout=20,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        raise AssetFailure(f"无法获取 Protenix 版本：{exc}") from exc
    text = completed.stdout.strip() or completed.stderr.strip()
    match = re.search(r"\b(\d+\.\d+\.\d+)\b", text)
    if completed.returncode != 0 or match is None:
        raise AssetFailure("Protenix CLI 未返回可校验的 x.y.z 版本。")
    return match.group(1)


def build_manifest(
    root: Path,
    model_name: str,
    protenix_version: str,
    *,
    include_template_caches: bool = False,
) -> dict[str, Any]:
    if not MODEL_RE.fullmatch(model_name):
        raise AssetFailure("model_name 不合法。")
    if root.is_symlink() or not root.is_dir():
        raise AssetFailure("protenix_root_dir 必须是已存在的非 symlink 目录。")
    root = root.resolve()
    checkpoint = checked_file(
        root,
        f"checkpoint/{model_name}.pt",
        "model_checkpoint",
    )
    caches = dict(BASE_CACHES)
    if include_template_caches:
        caches.update(TEMPLATE_CACHES)
    return {
        "schema_version": SCHEMA_VERSION,
        "model_name": model_name,
        "protenix_version": protenix_version,
        "protenix_root_dir": str(root),
        "checkpoint": checkpoint,
        "required_caches": [
            checked_file(root, relative, name)
            for name, relative in sorted(caches.items())
        ],
    }


def write_manifest(path: Path, document: dict[str, Any]) -> None:
    if path.exists() and path.is_symlink():
        raise AssetFailure("输出文件不得是 symlink。")
    temporary = path.with_name(f".{path.name}.tmp")
    temporary.write_text(
        json.dumps(document, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    temporary.replace(path)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="从已安装的 Protenix 资产构建固定 SHA-256 清单。"
    )
    parser.add_argument("--protenix-root-dir", required=True, type=Path)
    parser.add_argument("--model-name", required=True)
    parser.add_argument("--protenix-command", default="protenix")
    parser.add_argument("--protenix-version")
    parser.add_argument("--include-template-caches", action="store_true")
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        version = args.protenix_version or detect_version(args.protenix_command)
        if not re.fullmatch(r"\d+\.\d+\.\d+", version):
            raise AssetFailure("protenix_version 必须是 x.y.z。")
        document = build_manifest(
            args.protenix_root_dir,
            args.model_name,
            version,
            include_template_caches=args.include_template_caches,
        )
        args.output.parent.mkdir(parents=True, exist_ok=True)
        write_manifest(args.output, document)
    except (AssetFailure, OSError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 3
    print(f"wrote verified runtime assets manifest: {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
