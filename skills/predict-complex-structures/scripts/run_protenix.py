#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
import re
import shutil
import signal
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


SCHEMA_VERSION = "1.0.0"
WORKFLOW = "predict-complex-structures"
SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
BASE_CACHE_NAMES = {
    "ccd_components_file",
    "ccd_components_rdkit_mol_file",
    "obsolete_release_data_csv",
    "pdb_cluster_file",
}
TEMPLATE_CACHE_NAMES = {
    "obsolete_pdbs_path",
    "release_dates_path",
}
SECRET_RE = re.compile(
    r"(?i)(?:[\"']?(?:api[_ -]?key|access[_ -]?key|authorization|cookie|"
    r"credential|password|private[_ -]?key|secret|session[_ -]?key|token|"
    r"x-agent-plan-key)[\"']?\s*[:=]\s*[\"']?[^\"'\s,}]{4,}|"
    r"bearer\s+[A-Za-z0-9._~-]{16,}|"
    r"ark-[A-Za-z0-9-]{30,}|"
    r"sk-[A-Za-z0-9_-]{20,}|"
    r"AKIA[0-9A-Z]{16}|"
    r"github_pat_[A-Za-z0-9_]{20,}|"
    r"gh[pousr]_[A-Za-z0-9]{20,}|"
    r"eyJ[A-Za-z0-9_-]{20,}\.[A-Za-z0-9_-]{20,}\.[A-Za-z0-9_-]{20,})"
)
SENSITIVE_ENV_KEY_RE = re.compile(
    r"(?i)(?:(?:api|access|private|session)[_-]?key|authorization|cookie|"
    r"credential|password|secret|token)"
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
    "elapsed_seconds",
}
FINGERPRINT_KEYS = {
    "input_fingerprint",
    "evidence_fingerprint",
    "ligand_handoff_fingerprint",
    "run_fingerprint",
    "result_fingerprint",
}
_CANCEL_REQUESTED = False
JSON_LIMIT_BYTES = 64 * 1024 * 1024


class RunFailure(RuntimeError):
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
    payload = canonical_json(_without_nondeterministic_fields(document))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def value_sha256(value: Any) -> str:
    return hashlib.sha256(canonical_json(value).encode("utf-8")).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_json(path: Path, label: str) -> dict[str, Any]:
    if path.is_symlink() or not path.is_file():
        raise RunFailure(f"{label} 缺失或不是普通文件。")
    if path.stat().st_size > JSON_LIMIT_BYTES:
        raise RunFailure(f"{label} 超过 {JSON_LIMIT_BYTES} bytes。")
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (UnicodeError, json.JSONDecodeError) as exc:
        raise RunFailure(f"{label} 不是有效 UTF-8 JSON：{exc}") from exc
    if not isinstance(value, dict):
        raise RunFailure(f"{label} 顶层必须是 object。")
    return value


def read_json_value(path: Path, label: str) -> Any:
    if path.is_symlink() or not path.is_file():
        raise RunFailure(f"{label} 缺失或不是普通文件。")
    if path.stat().st_size > JSON_LIMIT_BYTES:
        raise RunFailure(f"{label} 超过 {JSON_LIMIT_BYTES} bytes。")
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (UnicodeError, json.JSONDecodeError) as exc:
        raise RunFailure(f"{label} 不是有效 UTF-8 JSON：{exc}") from exc


def verify_document_fingerprint(
    document: dict[str, Any],
    key: str,
    label: str,
) -> None:
    expected = document.get(key)
    if not isinstance(expected, str) or expected != document_fingerprint(document):
        raise RunFailure(f"{label}.{key} 缺失或与内容不一致。")


def validate_prepared_package(
    result_dir: Path,
    run_manifest: dict[str, Any],
) -> None:
    verify_document_fingerprint(
        run_manifest,
        "run_fingerprint",
        "run_manifest.json",
    )
    input_manifest = read_json(result_dir / "input_manifest.json", "input_manifest.json")
    if input_manifest.get("schema_version") != SCHEMA_VERSION:
        raise RunFailure(f"input_manifest.schema_version 必须是 {SCHEMA_VERSION}。")
    if input_manifest.get("workflow") != WORKFLOW:
        raise RunFailure("input_manifest.workflow 不正确。")
    verify_document_fingerprint(
        input_manifest,
        "input_fingerprint",
        "input_manifest.json",
    )
    if run_manifest.get("input_fingerprint") != input_manifest.get(
        "input_fingerprint"
    ):
        raise RunFailure("run_manifest 与 input_manifest 未绑定。")
    if run_manifest.get("job_id") != input_manifest.get("job_id"):
        raise RunFailure("run_manifest.job_id 与 input_manifest 不一致。")
    if run_manifest.get("model") != input_manifest.get("model_request"):
        raise RunFailure("run_manifest.model 与 input_manifest 未绑定。")

    protenix_input = read_json_value(
        result_dir / "protenix_input.json",
        "protenix_input.json",
    )
    if value_sha256(protenix_input) != input_manifest.get(
        "protenix_input_sha256"
    ):
        raise RunFailure("protenix_input.json 与 input_manifest 未绑定。")
    if run_manifest.get("inference_allowed") and not (
        isinstance(protenix_input, list) and protenix_input
    ):
        raise RunFailure("可执行任务的 protenix_input.json 必须是非空 array。")
    if not run_manifest.get("inference_allowed") and protenix_input != []:
        raise RunFailure("上游阻断时 protenix_input.json 必须是空 array。")

    for file_name, fingerprint_key, binding_key in (
        (
            "evidence.json",
            "evidence_fingerprint",
            "bound_evidence_fingerprint",
        ),
        (
            "ligand_handoff.json",
            "ligand_handoff_fingerprint",
            "bound_ligand_handoff_fingerprint",
        ),
    ):
        document = read_json(result_dir / file_name, file_name)
        if document.get("workflow") != WORKFLOW:
            raise RunFailure(f"{file_name}.workflow 不正确。")
        verify_document_fingerprint(document, fingerprint_key, file_name)
        if document.get(fingerprint_key) != input_manifest.get(binding_key):
            raise RunFailure(f"{file_name} 与 input_manifest 未绑定。")


def require_fresh_run_directory(
    result_dir: Path,
    run_manifest: dict[str, Any],
) -> None:
    if run_manifest.get("status") != "not_run":
        raise RunFailure("Runner 只接受 status=not_run 的全新准备结果包。")
    if run_manifest.get("command") is not None:
        raise RunFailure("结果包已记录运行命令；请重新 prepare 到新目录。")
    if run_manifest.get("started_at_utc") is not None:
        raise RunFailure("结果包已记录启动时间；请重新 prepare 到新目录。")
    outputs = run_manifest.get("outputs")
    if not isinstance(outputs, dict) or any(outputs.get(key) for key in (
        "structure_files",
        "auxiliary_structure_files",
        "summary_confidence_files",
        "full_confidence_files",
    )):
        raise RunFailure("结果包已登记预测输出；请重新 prepare 到新目录。")
    stale_names = (
        "protenix_output",
        "protenix_stdout.txt",
        "protenix_stderr.txt",
        "validation_report.json",
        "report.html",
    )
    stale = [name for name in stale_names if (result_dir / name).exists()]
    if stale:
        raise RunFailure(
            f"结果目录含陈旧运行产物 {stale}；请重新 prepare 到新目录。"
        )


def parse_cli_version(value: str) -> str | None:
    match = re.search(r"\b(\d+\.\d+\.\d+)\b", value)
    return match.group(1) if match else None


def validate_runtime_assets(
    manifest_path: Path,
    model: dict[str, Any],
) -> tuple[dict[str, Any], Path]:
    document = read_json(manifest_path, "runtime_assets_manifest")
    if document.get("schema_version") != SCHEMA_VERSION:
        raise RunFailure(
            f"runtime_assets_manifest.schema_version 必须是 {SCHEMA_VERSION}。"
        )
    if document.get("model_name") != model.get("name"):
        raise RunFailure("runtime_assets_manifest.model_name 与任务不一致。")
    root_raw = document.get("protenix_root_dir")
    if not isinstance(root_raw, str) or not root_raw:
        raise RunFailure("runtime_assets_manifest.protenix_root_dir 必填。")
    root = Path(root_raw)
    if not root.is_absolute() or root.is_symlink() or not root.is_dir():
        raise RunFailure("protenix_root_dir 必须是已存在的非 symlink 绝对目录。")
    root = root.resolve()

    def verify_asset(item: Any, label: str) -> dict[str, Any]:
        if not isinstance(item, dict):
            raise RunFailure(f"{label} 必须是 object。")
        relative_raw = item.get("relative_path")
        expected_hash = str(item.get("sha256") or "").lower()
        if not isinstance(relative_raw, str) or not relative_raw:
            raise RunFailure(f"{label}.relative_path 必填。")
        relative = Path(relative_raw)
        if relative.is_absolute() or ".." in relative.parts:
            raise RunFailure(f"{label}.relative_path 必须位于 Protenix 根目录内。")
        if not SHA256_RE.fullmatch(expected_hash):
            raise RunFailure(f"{label}.sha256 必须是 64 位小写十六进制。")
        unresolved_path = root / relative
        if unresolved_path.is_symlink():
            raise RunFailure(f"{label} 不允许 symlink。")
        path = unresolved_path.resolve()
        try:
            path.relative_to(root)
        except ValueError as exc:
            raise RunFailure(f"{label} 路径逃逸 Protenix 根目录。") from exc
        if path.is_symlink() or not path.is_file() or path.stat().st_size <= 0:
            raise RunFailure(f"{label} 文件缺失、为空或为 symlink。")
        actual_hash = sha256_file(path)
        if actual_hash != expected_hash:
            raise RunFailure(f"{label} SHA-256 与资产清单不一致。")
        return {
            "name": item.get("name"),
            "relative_path": relative.as_posix(),
            "sha256": actual_hash,
            "size_bytes": path.stat().st_size,
        }

    checkpoint = verify_asset(document.get("checkpoint"), "checkpoint")
    if Path(checkpoint["relative_path"]).name != f"{model.get('name')}.pt":
        raise RunFailure("checkpoint 文件名必须与 model.name 对应。")
    raw_caches = document.get("required_caches")
    if not isinstance(raw_caches, list):
        raise RunFailure("runtime_assets_manifest.required_caches 必须是 array。")
    caches = [
        verify_asset(item, f"required_caches[{index}]")
        for index, item in enumerate(raw_caches)
    ]
    cache_names = {item.get("name") for item in caches}
    required_names = set(BASE_CACHE_NAMES)
    if model.get("template_strategy") != "none":
        required_names.update(TEMPLATE_CACHE_NAMES)
    missing = sorted(required_names - cache_names)
    if missing:
        raise RunFailure(f"runtime_assets_manifest 缺少运行缓存：{missing}。")
    summary = {
        "manifest_sha256": sha256_file(manifest_path),
        "model_name": model.get("name"),
        "protenix_version": document.get("protenix_version"),
        "checkpoint": checkpoint,
        "required_caches": sorted(
            caches,
            key=lambda item: str(item.get("name")),
        ),
        "downloads_allowed": False,
    }
    return summary, root


def write_json(path: Path, document: dict[str, Any]) -> None:
    text = json.dumps(document, ensure_ascii=False, indent=2) + "\n"
    if SECRET_RE.search(text):
        raise RunFailure(f"{path.name} 含疑似凭证。")
    if SENSITIVE_PATH_RE.search(text):
        raise RunFailure(f"{path.name} 含主机敏感路径。")
    temporary = path.with_name(f".{path.name}.tmp")
    temporary.write_text(text, encoding="utf-8")
    temporary.replace(path)


def write_text_safely(path: Path, text: str) -> None:
    if path.is_symlink():
        raise RunFailure(f"{path.name} 不允许 symlink。")
    temporary = path.with_name(f".{path.name}.tmp")
    temporary.write_text(text, encoding="utf-8")
    temporary.replace(path)


def sanitize_text(text: str, result_dir: Path) -> str:
    sanitized = text
    replacements = {
        str(result_dir.resolve()): "<RESULT_DIR>",
        str(result_dir.resolve().parent): "<WORK_DIR>",
        str(Path.home()): "<HOME>",
    }
    for original, replacement in sorted(
        replacements.items(),
        key=lambda item: len(item[0]),
        reverse=True,
    ):
        if original:
            sanitized = sanitized.replace(original, replacement)
    sanitized = SECRET_RE.sub("<REDACTED_CREDENTIAL>", sanitized)
    sanitized = SENSITIVE_PATH_RE.sub("<REDACTED_PATH>", sanitized)
    return sanitized


def sanitized_process_environment() -> tuple[dict[str, str], list[str]]:
    environment: dict[str, str] = {}
    removed: list[str] = []
    for key, value in os.environ.items():
        if SENSITIVE_ENV_KEY_RE.search(key):
            removed.append(key)
            continue
        environment[key] = value
    return environment, sorted(removed)


def load_validator() -> Any:
    path = Path(__file__).with_name("validate_output.py")
    spec = importlib.util.spec_from_file_location(
        "predict_complex_structures_validator",
        path,
    )
    if spec is None or spec.loader is None:
        raise RunFailure("无法加载独立 Validator。")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def capture_command(command: list[str], timeout: float = 20.0) -> dict[str, Any]:
    try:
        completed = subprocess.run(
            command,
            capture_output=True,
            text=True,
            check=False,
            timeout=timeout,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        return {
            "returncode": None,
            "stdout": "",
            "stderr": str(exc),
        }
    return {
        "returncode": completed.returncode,
        "stdout": completed.stdout,
        "stderr": completed.stderr,
    }


def detect_gpu(nvidia_smi_command: str = "nvidia-smi") -> dict[str, Any] | None:
    executable = shutil.which(nvidia_smi_command)
    if executable is None:
        return None
    result = capture_command(
        [
            executable,
            "--query-gpu=name,memory.total,driver_version",
            "--format=csv,noheader,nounits",
        ]
    )
    if result["returncode"] != 0:
        return None
    devices = []
    for line in result["stdout"].splitlines():
        values = [item.strip() for item in line.split(",")]
        if len(values) >= 3:
            devices.append(
                {
                    "name": values[0],
                    "memory_total_mb": values[1],
                    "driver_version": values[2],
                }
            )
    return {"devices": devices} if devices else None


def build_command(
    executable: str,
    manifest: dict[str, Any],
) -> list[str]:
    model = manifest.get("model")
    if not isinstance(model, dict):
        raise RunFailure("run_manifest.model 缺失。")
    seeds = model.get("seeds")
    if not isinstance(seeds, list) or not seeds:
        raise RunFailure("run_manifest.model.seeds 缺失。")
    command = [
        executable,
        "pred",
        "-i",
        "protenix_input.json",
        "-o",
        "protenix_output",
        "-s",
        ",".join(str(item) for item in seeds),
        "-n",
        str(model.get("name")),
        "--use_msa",
        "true" if model.get("msa_strategy") != "none" else "false",
        "--use_template",
        "true" if model.get("template_strategy") != "none" else "false",
        "--use_rna_msa",
        "true" if model.get("rna_msa_strategy") != "none" else "false",
        "--dtype",
        str(model.get("dtype", "bf16")),
        "--sample",
        str(model.get("samples", 5)),
        "--step",
        str(model.get("steps", 200)),
        "--cycle",
        str(model.get("cycles", 10)),
    ]
    return command


def terminate_process_group(process: subprocess.Popen[str]) -> None:
    if process.poll() is not None:
        return
    try:
        if os.name == "posix":
            os.killpg(process.pid, signal.SIGTERM)
        else:
            process.terminate()
        process.wait(timeout=5)
        return
    except (OSError, subprocess.TimeoutExpired):
        pass
    try:
        if os.name == "posix":
            os.killpg(process.pid, signal.SIGKILL)
        else:
            process.kill()
        process.wait(timeout=5)
    except (OSError, subprocess.TimeoutExpired):
        process.kill()


def execute_process(
    command: list[str],
    cwd: Path,
    timeout_seconds: float,
    cancel_file: Path | None = None,
    environment: dict[str, str] | None = None,
) -> dict[str, Any]:
    global _CANCEL_REQUESTED
    started = time.monotonic()
    try:
        process = subprocess.Popen(
            command,
            cwd=cwd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            start_new_session=(os.name == "posix"),
            env=environment,
        )
    except OSError as exc:
        return {
            "exit_code": None,
            "stdout": "",
            "stderr": str(exc),
            "timed_out": False,
            "canceled": False,
            "elapsed_seconds": time.monotonic() - started,
        }
    timed_out = False
    canceled = False
    while True:
        try:
            stdout, stderr = process.communicate(timeout=0.25)
            break
        except subprocess.TimeoutExpired:
            if _CANCEL_REQUESTED or (
                cancel_file is not None and cancel_file.exists()
            ):
                canceled = True
                terminate_process_group(process)
                stdout, stderr = process.communicate()
                break
            if time.monotonic() - started >= timeout_seconds:
                timed_out = True
                terminate_process_group(process)
                stdout, stderr = process.communicate()
                break
    return {
        "exit_code": process.returncode,
        "stdout": stdout,
        "stderr": stderr,
        "timed_out": timed_out,
        "canceled": canceled,
        "elapsed_seconds": time.monotonic() - started,
    }


def discover_outputs(result_dir: Path) -> dict[str, list[str]]:
    output_root = result_dir / "protenix_output"
    structures: list[str] = []
    auxiliary_structures: list[str] = []
    summaries: list[str] = []
    full: list[str] = []
    if output_root.is_dir():
        for path in sorted(output_root.rglob("*")):
            if path.is_symlink() or not path.is_file():
                continue
            relative = path.relative_to(result_dir).as_posix()
            if path.suffix.lower() in {".cif", ".mmcif"}:
                if path.stem.endswith("_wounresol"):
                    auxiliary_structures.append(relative)
                else:
                    structures.append(relative)
            elif (
                path.suffix.lower() == ".json"
                and "_summary_confidence_sample_" in path.name
            ):
                summaries.append(relative)
            elif (
                path.suffix.lower() == ".json"
                and "_full_data_sample_" in path.name
            ):
                full.append(relative)
    return {
        "structure_files": structures,
        "auxiliary_structure_files": auxiliary_structures,
        "summary_confidence_files": summaries,
        "full_confidence_files": full,
    }


def update_manifest(
    path: Path,
    manifest: dict[str, Any],
) -> None:
    manifest["run_fingerprint"] = document_fingerprint(manifest)
    write_json(path, manifest)


def write_logs(
    result_dir: Path,
    stdout: str,
    stderr: str,
) -> None:
    write_text_safely(
        result_dir / "protenix_stdout.txt",
        sanitize_text(stdout, result_dir),
    )
    write_text_safely(
        result_dir / "protenix_stderr.txt",
        sanitize_text(stderr, result_dir),
    )


def validate_and_write(result_dir: Path, generated_at_utc: str | None) -> dict[str, Any]:
    validator = load_validator()
    return validator.validate_result_dir(
        result_dir,
        generated_at_utc=generated_at_utc,
        write_artifacts=True,
    )


def stop_without_running(
    result_dir: Path,
    manifest: dict[str, Any],
    manifest_path: Path,
    *,
    failure_code: str,
    message: str,
    generated_at_utc: str | None,
    exit_code: int = 2,
) -> tuple[dict[str, Any], int]:
    manifest.update(
        {
            "status": "not_run",
            "finished_at_utc": generated_at_utc or now_utc(),
            "exit_code": None,
            "failure": {"code": failure_code, "message": message},
            "outputs": {
                "structure_files": [],
                "auxiliary_structure_files": [],
                "summary_confidence_files": [],
                "full_confidence_files": [],
            },
        }
    )
    update_manifest(manifest_path, manifest)
    write_logs(result_dir, "", message + "\n")
    report = validate_and_write(result_dir, generated_at_utc)
    return report, exit_code


def run_job(
    result_dir: Path,
    *,
    protenix_command: str = "protenix",
    nvidia_smi_command: str = "nvidia-smi",
    timeout_seconds: float = 3600,
    cancel_file: Path | None = None,
    runtime_assets_manifest: Path | None = None,
    allow_download: bool = False,
    generated_at_utc: str | None = None,
) -> tuple[dict[str, Any], int]:
    if result_dir.is_symlink():
        raise RunFailure("result_dir 不允许 symlink。")
    result_dir = result_dir.resolve()
    manifest_path = result_dir / "run_manifest.json"
    manifest = read_json(manifest_path, "run_manifest.json")
    if manifest.get("schema_version") != SCHEMA_VERSION:
        raise RunFailure(f"schema_version 必须是 {SCHEMA_VERSION}。")
    if manifest.get("workflow") != WORKFLOW:
        raise RunFailure("run_manifest.workflow 不正确。")
    validate_prepared_package(result_dir, manifest)
    if not manifest.get("inference_allowed"):
        return stop_without_running(
            result_dir,
            manifest,
            manifest_path,
            failure_code="UPSTREAM_BLOCKED",
            message="上游身份或结构状态阻止 Protenix 推理。",
            generated_at_utc=generated_at_utc,
        )
    require_fresh_run_directory(result_dir, manifest)
    if timeout_seconds <= 0:
        raise RunFailure("timeout_seconds 必须大于 0。")

    executable = shutil.which(protenix_command)
    if executable is None:
        return stop_without_running(
            result_dir,
            manifest,
            manifest_path,
            failure_code="PROTENIX_NOT_INSTALLED",
            message="未找到官方 protenix CLI；未执行结构预测。",
            generated_at_utc=generated_at_utc,
        )
    gpu = detect_gpu(nvidia_smi_command)
    if gpu is None:
        return stop_without_running(
            result_dir,
            manifest,
            manifest_path,
            failure_code="GPU_UNAVAILABLE",
            message="未发现可用 NVIDIA GPU；未执行 Protenix。",
            generated_at_utc=generated_at_utc,
        )

    version_result = capture_command([executable, "--version"])
    version = (
        version_result["stdout"].strip()
        or version_result["stderr"].strip()
        or "not_exposed_by_cli"
    )
    actual_version = parse_cli_version(version)
    model = manifest.get("model")
    if not isinstance(model, dict):
        raise RunFailure("run_manifest.model 缺失。")
    requested_version = str(model.get("version") or "")
    if (
        requested_version not in {"", "not_exposed_by_cli", "unverified"}
        and actual_version != requested_version
    ):
        return stop_without_running(
            result_dir,
            manifest,
            manifest_path,
            failure_code="RUNTIME_VERSION_MISMATCH",
            message=(
                f"Protenix 版本不一致：requested={requested_version}, "
                f"actual={actual_version or 'unknown'}。"
            ),
            generated_at_utc=generated_at_utc,
        )

    asset_summary: dict[str, Any]
    process_environment, removed_environment_keys = (
        sanitized_process_environment()
    )
    if allow_download:
        asset_summary = {
            "mode": "download_explicitly_allowed",
            "downloads_allowed": True,
        }
    else:
        if runtime_assets_manifest is None:
            return stop_without_running(
                result_dir,
                manifest,
                manifest_path,
                failure_code="ASSETS_MANIFEST_REQUIRED",
                message=(
                    "默认禁止隐式下载；必须提供已校验的 runtime assets manifest，"
                    "或经用户确认显式使用 --allow-download。"
                ),
                generated_at_utc=generated_at_utc,
            )
        try:
            asset_summary, protenix_root = validate_runtime_assets(
                runtime_assets_manifest,
                model,
            )
        except RunFailure as exc:
            return stop_without_running(
                result_dir,
                manifest,
                manifest_path,
                failure_code="RUNTIME_ASSETS_INVALID",
                message=str(exc),
                generated_at_utc=generated_at_utc,
            )
        expected_version = asset_summary.get("protenix_version")
        if expected_version and expected_version != actual_version:
            return stop_without_running(
                result_dir,
                manifest,
                manifest_path,
                failure_code="RUNTIME_VERSION_MISMATCH",
                message=(
                    f"资产清单 Protenix 版本与 CLI 不一致："
                    f"manifest={expected_version}, actual={actual_version or 'unknown'}。"
                ),
                generated_at_utc=generated_at_utc,
            )
        process_environment["PROTENIX_ROOT_DIR"] = str(protenix_root)

    command = build_command(executable, manifest)
    display_command = [Path(command[0]).name, *command[1:]]
    started_at = generated_at_utc or now_utc()
    manifest.update(
        {
            "status": "partial",
            "runtime": {
                "backend": "protenix_cli",
                "environment": "independent_gpu_runtime",
                "protenix_version": actual_version or sanitize_text(
                    version,
                    result_dir,
                ),
                "device": gpu,
                "assets": asset_summary,
                "removed_sensitive_environment_key_count": len(
                    removed_environment_keys
                ),
            },
            "command": display_command,
            "started_at_utc": started_at,
            "finished_at_utc": None,
            "elapsed_seconds": None,
            "exit_code": None,
            "failure": None,
        }
    )
    update_manifest(manifest_path, manifest)

    process_result = execute_process(
        command,
        result_dir,
        timeout_seconds,
        cancel_file,
        process_environment,
    )
    write_logs(
        result_dir,
        process_result["stdout"],
        process_result["stderr"],
    )
    manifest["finished_at_utc"] = generated_at_utc or now_utc()
    manifest["elapsed_seconds"] = round(
        float(process_result["elapsed_seconds"]),
        6,
    )
    manifest["exit_code"] = process_result["exit_code"]
    manifest["outputs"] = discover_outputs(result_dir)

    if process_result["canceled"]:
        manifest["status"] = "error"
        manifest["failure"] = {
            "code": "CANCELED",
            "message": "用户取消了 Protenix 任务。",
        }
        update_manifest(manifest_path, manifest)
        report = validate_and_write(result_dir, generated_at_utc)
        return report, 130
    if process_result["timed_out"]:
        manifest["status"] = "error"
        manifest["failure"] = {
            "code": "TIMEOUT",
            "message": f"Protenix 超过 {timeout_seconds} 秒。",
        }
        update_manifest(manifest_path, manifest)
        report = validate_and_write(result_dir, generated_at_utc)
        return report, 1
    if process_result["exit_code"] != 0:
        manifest["status"] = "error"
        manifest["failure"] = {
            "code": "NONZERO_EXIT",
            "message": (
                f"Protenix 返回非零退出码 {process_result['exit_code']}。"
            ),
        }
        update_manifest(manifest_path, manifest)
        report = validate_and_write(result_dir, generated_at_utc)
        return report, 1

    update_manifest(manifest_path, manifest)
    provisional = validate_and_write(result_dir, generated_at_utc)
    if not provisional.get("completion_ready"):
        manifest["status"] = "error"
        manifest["failure"] = {
            "code": "OUTPUT_INVALID",
            "message": "Protenix 输出缺失或未通过独立 Validator。",
        }
        update_manifest(manifest_path, manifest)
        report = validate_and_write(result_dir, generated_at_utc)
        return report, 1

    manifest["status"] = "completed"
    manifest["review_disposition"] = provisional["review_disposition"]
    manifest["failure"] = None
    update_manifest(manifest_path, manifest)
    final_report = validate_and_write(result_dir, generated_at_utc)
    if not final_report.get("valid"):
        manifest["status"] = "error"
        manifest["review_disposition"] = "rejected"
        manifest["failure"] = {
            "code": "OUTPUT_INVALID",
            "message": "completed 状态未通过最终 Validator。",
        }
        update_manifest(manifest_path, manifest)
        final_report = validate_and_write(result_dir, generated_at_utc)
        return final_report, 1
    return final_report, 0


def _handle_cancel_signal(signum: int, _frame: Any) -> None:
    global _CANCEL_REQUESTED
    _CANCEL_REQUESTED = True


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="在独立 GPU Runtime 中执行官方 Protenix CLI。"
    )
    parser.add_argument("--result-dir", required=True, type=Path)
    parser.add_argument("--protenix-command", default="protenix")
    parser.add_argument("--nvidia-smi-command", default="nvidia-smi")
    parser.add_argument("--timeout-seconds", type=float, default=3600)
    parser.add_argument("--cancel-file", type=Path)
    parser.add_argument("--runtime-assets-manifest", type=Path)
    parser.add_argument(
        "--allow-download",
        action="store_true",
        help="经用户明确确认后，允许 Protenix 下载缺失权重和缓存。",
    )
    parser.add_argument("--generated-at")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    signal.signal(signal.SIGTERM, _handle_cancel_signal)
    signal.signal(signal.SIGINT, _handle_cancel_signal)
    try:
        report, exit_code = run_job(
            args.result_dir,
            protenix_command=args.protenix_command,
            nvidia_smi_command=args.nvidia_smi_command,
            timeout_seconds=args.timeout_seconds,
            cancel_file=args.cancel_file,
            runtime_assets_manifest=args.runtime_assets_manifest,
            allow_download=args.allow_download,
            generated_at_utc=args.generated_at,
        )
    except (OSError, RunFailure, RuntimeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 3
    print(json.dumps(report, ensure_ascii=False, indent=2))
    return exit_code


if __name__ == "__main__":
    raise SystemExit(main())
