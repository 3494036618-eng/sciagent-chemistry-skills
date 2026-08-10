#!/usr/bin/env python3
"""保守解析化学名称、结构和数据库标识符，并保留来源与不确定性。"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import re
import socket
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable, Iterable, Optional, Sequence


SCHEMA_VERSION = "1.0.0"
WORKFLOW = "chemical-identity-resolution"
USER_AGENT = "resolve-chemical-identities/1.0"
DEFAULT_TIMEOUT = 20
DEFAULT_RETRIES = 1
DEFAULT_SOURCES = ("opsin", "pubchem", "chembl", "unichem")
SUPPORTED_SOURCES = frozenset(DEFAULT_SOURCES)
INPUT_TYPES = frozenset(
    {
        "auto",
        "name",
        "smiles",
        "inchi",
        "inchikey",
        "pubchem_cid",
        "chembl_id",
        "cas_rn",
    }
)
INPUT_TYPE_ALIASES = {
    "cid": "pubchem_cid",
    "pubchem-cid": "pubchem_cid",
    "chembl-id": "chembl_id",
    "cas": "cas_rn",
    "cas-rn": "cas_rn",
}
INCHIKEY_RE = re.compile(r"^[A-Z]{14}-[A-Z]{10}-[A-Z]$")
CHEMBL_RE = re.compile(r"^CHEMBL\d+$", re.IGNORECASE)
CAS_RE = re.compile(r"^(?P<body>\d{2,7})-(?P<middle>\d{2})-(?P<check>\d)$")
SECRET_RE = re.compile(
    r"(?i)(authorization\s*:|bearer\s+[A-Za-z0-9._~+/=-]{12,}|"
    r"api[_ -]?key\s*[:=]|cookie\s*:|ark-[A-Za-z0-9_-]{16,})"
)
PUBCHEM_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
CHEMBL_BASE = "https://www.ebi.ac.uk/chembl/api/data"
OPSIN_BASE = "https://www.ebi.ac.uk/opsin/ws"
UNICHEM_BASE = "https://www.ebi.ac.uk/unichem/api/v1"
PUBCHEM_PROPERTIES = (
    "Title,IUPACName,MolecularFormula,CanonicalSMILES,IsomericSMILES,"
    "ConnectivitySMILES,InChI,InChIKey"
)
TEMPORAL_KEYS = frozenset(
    {
        "generated_at_utc",
        "retrieved_at_utc",
        "requested_at_utc",
        "confirmed_at_utc",
    }
)


class DependencyFailure(RuntimeError):
    """固定依赖不存在或不可加载。"""


class InputFailure(ValueError):
    """请求文件或命令行输入无法加载。"""


def now_utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def canonical_json(value: Any) -> str:
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"))


def sha256_json(value: Any) -> str:
    return hashlib.sha256(canonical_json(value).encode("utf-8")).hexdigest()


def _without_temporal_fields(value: Any) -> Any:
    if isinstance(value, dict):
        return {
            key: _without_temporal_fields(item)
            for key, item in value.items()
            if key not in TEMPORAL_KEYS and key != "result_fingerprint"
        }
    if isinstance(value, list):
        return [_without_temporal_fields(item) for item in value]
    return value


def output_fingerprint(document: dict[str, Any]) -> str:
    return sha256_json(_without_temporal_fields(document))


def load_toolkit() -> dict[str, Any]:
    try:
        import rdkit
        from rdkit import Chem, RDLogger
        from rdkit.Chem import inchi, rdMolDescriptors
    except ImportError as error:
        raise DependencyFailure(
            "缺少 rdkit==2025.9.2；请在隔离环境安装 scripts/requirements.txt"
        ) from error

    RDLogger.DisableLog("rdApp.error")
    return {
        "rdkit": rdkit,
        "Chem": Chem,
        "inchi": inchi,
        "rdMolDescriptors": rdMolDescriptors,
    }


def toolkit_versions(toolkit: dict[str, Any]) -> dict[str, Any]:
    return {
        "resolver": "1.0.0",
        "rdkit": toolkit["rdkit"].__version__,
        "inchi_provider": {
            "name": "rdkit.Chem.inchi",
            "embedded_inchi_version": "not_exposed_by_rdkit_python_api",
        },
        "opsin": {
            "runtime": "official_web_api_when_enabled",
            "service_version": "not_exposed_by_api",
            "reviewed_release": "2.9.0",
        },
        "pubchem": {"api": "PUG REST", "service_version": "not_exposed_by_api"},
        "chembl": {"api": "ChEMBL Data Web Services", "database_version": "runtime"},
        "unichem": {"api": "POST API v1 (2.0 documentation)"},
    }


def normalize_input_type(value: str) -> str:
    normalized = value.strip().lower().replace(" ", "_")
    normalized = INPUT_TYPE_ALIASES.get(normalized, normalized)
    if normalized not in INPUT_TYPES:
        raise InputFailure(f"不支持的 input_type：{value}")
    return normalized


def valid_cas_check_digit(value: str) -> bool:
    match = CAS_RE.fullmatch(value)
    if not match:
        return False
    digits = match.group("body") + match.group("middle")
    total = sum(int(digit) * weight for weight, digit in enumerate(reversed(digits), 1))
    return total % 10 == int(match.group("check"))


def parse_structure(
    value: str, structure_type: str, toolkit: dict[str, Any]
) -> tuple[Optional[Any], Optional[str]]:
    Chem = toolkit["Chem"]
    try:
        if structure_type == "smiles":
            molecule = Chem.MolFromSmiles(value, sanitize=False)
        elif structure_type == "inchi":
            molecule = Chem.MolFromInchi(value, sanitize=False, removeHs=False)
        else:
            raise ValueError(f"不支持的结构类型：{structure_type}")
        if molecule is None:
            return None, "RDKit 未生成分子对象"
        Chem.SanitizeMol(molecule)
        return molecule, None
    except Exception as error:  # RDKit raises several C++-backed exception types.
        return None, str(error)


def structure_identifiers(molecule: Any, toolkit: dict[str, Any]) -> dict[str, Any]:
    Chem = toolkit["Chem"]
    inchi = toolkit["inchi"]
    rdMolDescriptors = toolkit["rdMolDescriptors"]
    canonical_smiles = Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=True)
    standard_inchi = inchi.MolToInchi(molecule)
    standard_inchikey = inchi.InchiToInchiKey(standard_inchi) if standard_inchi else None
    chiral_centers = Chem.FindMolChiralCenters(
        molecule, includeUnassigned=True, includeCIP=True
    )
    unassigned_stereo = [
        {"atom_index": atom_index, "assignment": assignment}
        for atom_index, assignment in chiral_centers
        if assignment == "?"
    ]
    return {
        "canonical_smiles": canonical_smiles,
        "inchi": standard_inchi or None,
        "inchikey": standard_inchikey or None,
        "connectivity_block": (
            standard_inchikey[:14] if standard_inchikey else None
        ),
        "molecular_formula": rdMolDescriptors.CalcMolFormula(molecule),
        "component_count": len(Chem.GetMolFrags(molecule)),
        "unassigned_stereo": unassigned_stereo,
    }


def looks_like_failed_smiles(value: str) -> bool:
    if any(character.isspace() for character in value):
        return False
    if any(token in value for token in ("[", "]", "(", ")", "=", "#", "@", "\\", "/")):
        return True
    return bool(re.search(r"[A-Za-z]\d|\d[A-Za-z]", value))


def detect_input_type(
    query: str, toolkit: dict[str, Any]
) -> tuple[str, list[dict[str, Any]]]:
    stripped = query.strip()
    findings: list[dict[str, Any]] = []
    if CHEMBL_RE.fullmatch(stripped):
        return "chembl_id", findings
    if INCHIKEY_RE.fullmatch(stripped.upper()):
        return "inchikey", findings
    if stripped.startswith("InChI="):
        return "inchi", findings
    if CAS_RE.fullmatch(stripped):
        return "cas_rn", findings
    if stripped.isdigit():
        findings.append(
            {
                "code": "E-AMBIGUOUS-NUMERIC-ID",
                "severity": "error",
                "message": (
                    "纯数字不能安全地区分 PubChem CID、内部编号或名称；"
                    "请显式指定 input_type=pubchem_cid。"
                ),
            }
        )
        return "ambiguous_numeric", findings

    molecule, _ = parse_structure(stripped, "smiles", toolkit)
    if molecule is not None:
        return "smiles", findings
    if looks_like_failed_smiles(stripped):
        return "smiles", findings
    return "name", findings


def local_source_record(
    request_id: str,
    query: str,
    input_type: str,
    toolkit: dict[str, Any],
) -> tuple[Optional[dict[str, Any]], list[dict[str, Any]]]:
    molecule, error = parse_structure(query, input_type, toolkit)
    if molecule is None:
        return None, [
            {
                "code": "E-INVALID-STRUCTURE",
                "severity": "error",
                "message": error or "结构无法解析",
            }
        ]

    identifiers = structure_identifiers(molecule, toolkit)
    findings: list[dict[str, Any]] = []
    if identifiers["component_count"] > 1:
        findings.append(
            {
                "code": "R-MULTICOMPONENT-INPUT",
                "severity": "review",
                "message": "输入包含多个结构组分；不得自动挑选单一主体。",
            }
        )
    if identifiers["unassigned_stereo"]:
        findings.append(
            {
                "code": "R-UNSPECIFIED-STEREO",
                "severity": "review",
                "message": "结构存在未指定立体中心；不得自动补成立体化学。",
            }
        )
    return (
        {
            "source": "local_input",
            "source_family": "local_input",
            "source_record_id": request_id,
            "match_method": f"parsed_{input_type}",
            "title": None,
            "names": [],
            "structure": query if input_type == "smiles" else identifiers["canonical_smiles"],
            "inchi": identifiers["inchi"],
            "inchikey": identifiers["inchikey"],
            "molecular_formula": identifiers["molecular_formula"],
            "source_url": None,
            "raw_record": {
                "input_type": input_type,
                "input_value": query,
            },
            "record_findings": findings,
        },
        findings,
    )


def validate_request(
    item: dict[str, Any], toolkit: dict[str, Any]
) -> dict[str, Any]:
    request_id = str(item.get("id") or "query-1").strip()
    raw_query = item.get("query")
    query = raw_query if isinstance(raw_query, str) else ""
    requested_type = normalize_input_type(str(item.get("input_type") or "auto"))
    context = item.get("context")
    expected_form = item.get("expected_form")
    findings: list[dict[str, Any]] = []

    if not request_id:
        request_id = "query-1"
    if not query.strip():
        findings.append(
            {
                "code": "E-EMPTY-QUERY",
                "severity": "error",
                "message": "query 不能为空。",
            }
        )
        return {
            "id": request_id,
            "query": query,
            "requested_input_type": requested_type,
            "detected_input_type": "unknown",
            "context": context,
            "expected_form": expected_form,
            "input_status": "invalid_input",
            "findings": findings,
            "local_record": None,
        }

    stripped = query.strip()
    if len(stripped) > 1000:
        findings.append(
            {
                "code": "E-QUERY-TOO-LONG",
                "severity": "error",
                "message": "单条 query 超过 1000 个字符；请确认没有误传整段文档。",
            }
        )
        detected_type = "unknown"
    elif requested_type == "auto":
        detected_type, detected_findings = detect_input_type(stripped, toolkit)
        findings.extend(detected_findings)
    else:
        detected_type = requested_type

    local_record = None
    if detected_type == "ambiguous_numeric":
        pass
    elif detected_type == "inchikey":
        normalized_key = stripped.upper()
        if not INCHIKEY_RE.fullmatch(normalized_key):
            findings.append(
                {
                    "code": "E-INVALID-INCHIKEY",
                    "severity": "error",
                    "message": "InChIKey 必须符合 14-10-1 大写字母格式。",
                }
            )
    elif detected_type == "chembl_id":
        if not CHEMBL_RE.fullmatch(stripped):
            findings.append(
                {
                    "code": "E-INVALID-CHEMBL-ID",
                    "severity": "error",
                    "message": "ChEMBL ID 必须符合 CHEMBL 加数字的格式。",
                }
            )
    elif detected_type == "pubchem_cid":
        if not stripped.isdigit() or int(stripped) <= 0:
            findings.append(
                {
                    "code": "E-INVALID-PUBCHEM-CID",
                    "severity": "error",
                    "message": "PubChem CID 必须是正整数。",
                }
            )
    elif detected_type == "cas_rn":
        if not valid_cas_check_digit(stripped):
            findings.append(
                {
                    "code": "E-INVALID-CAS-CHECK-DIGIT",
                    "severity": "error",
                    "message": (
                        "CAS RN 格式或校验位无效；本检查只验证字符和校验位，"
                        "不证明该编号由 CAS 正式登记。"
                    ),
                }
            )
    elif detected_type in {"smiles", "inchi"}:
        local_record, local_findings = local_source_record(
            request_id, stripped, detected_type, toolkit
        )
        findings.extend(local_findings)
    elif detected_type == "name":
        if len(stripped) > 500:
            findings.append(
                {
                    "code": "E-NAME-TOO-LONG",
                    "severity": "error",
                    "message": "化学名称超过 500 个字符；请检查输入。",
                }
            )

    invalid = any(finding["severity"] == "error" for finding in findings)
    return {
        "id": request_id,
        "query": query,
        "normalized_query": (
            stripped.upper()
            if detected_type in {"inchikey", "chembl_id"}
            else stripped
        ),
        "requested_input_type": requested_type,
        "detected_input_type": detected_type,
        "context": context,
        "expected_form": expected_form,
        "input_status": "invalid_input" if invalid else "valid",
        "findings": findings,
        "local_record": local_record,
    }


class HttpTransport:
    def __init__(
        self,
        timeout: int = DEFAULT_TIMEOUT,
        retries: int = DEFAULT_RETRIES,
        clock: Callable[[], str] = now_utc,
        sleep: Callable[[float], None] = time.sleep,
    ):
        self.timeout = timeout
        self.retries = retries
        self.clock = clock
        self.sleep = sleep

    def request_json(
        self,
        fixture_key: str,
        method: str,
        url: str,
        body: Optional[dict[str, Any]] = None,
        body_format: str = "json",
    ) -> dict[str, Any]:
        del fixture_key
        requested_at = self.clock()
        encoded_body: Optional[bytes] = None
        headers = {
            "Accept": "application/json",
            "User-Agent": USER_AGENT,
        }
        if body is not None:
            if body_format == "form":
                encoded_body = urllib.parse.urlencode(body).encode("utf-8")
                headers["Content-Type"] = "application/x-www-form-urlencoded"
            else:
                encoded_body = canonical_json(body).encode("utf-8")
                headers["Content-Type"] = "application/json"
        request_body_sha256 = (
            hashlib.sha256(encoded_body).hexdigest() if encoded_body else None
        )

        for attempt in range(self.retries + 1):
            request = urllib.request.Request(
                url=url,
                data=encoded_body,
                headers=headers,
                method=method,
            )
            try:
                with urllib.request.urlopen(request, timeout=self.timeout) as response:
                    raw = response.read()
                    http_status = int(response.status)
                try:
                    payload = json.loads(raw.decode("utf-8"))
                except (UnicodeDecodeError, json.JSONDecodeError) as error:
                    return self._result(
                        requested_at,
                        method,
                        url,
                        "source_error",
                        http_status,
                        "invalid_json",
                        str(error),
                        None,
                        raw,
                        request_body_sha256,
                    )
                return self._result(
                    requested_at,
                    method,
                    url,
                    "success",
                    http_status,
                    None,
                    None,
                    payload,
                    raw,
                    request_body_sha256,
                )
            except urllib.error.HTTPError as error:
                raw = error.read()
                http_status = int(error.code)
                payload = None
                try:
                    payload = json.loads(raw.decode("utf-8"))
                except (UnicodeDecodeError, json.JSONDecodeError):
                    pass
                if http_status == 404:
                    return self._result(
                        requested_at,
                        method,
                        url,
                        "not_found",
                        http_status,
                        "not_found",
                        _payload_message(payload) or "HTTP 404",
                        payload,
                        raw,
                        request_body_sha256,
                    )
                retryable = http_status in {429, 500, 502, 503, 504}
                if retryable and attempt < self.retries:
                    self.sleep(min(2**attempt, 2))
                    continue
                error_kind = (
                    "rate_limited"
                    if http_status == 429
                    else "service_error"
                    if retryable
                    else "http_error"
                )
                return self._result(
                    requested_at,
                    method,
                    url,
                    "source_error",
                    http_status,
                    error_kind,
                    _payload_message(payload) or f"HTTP {http_status}",
                    payload,
                    raw,
                    request_body_sha256,
                )
            except (urllib.error.URLError, socket.timeout, TimeoutError) as error:
                if attempt < self.retries:
                    self.sleep(min(2**attempt, 2))
                    continue
                reason = getattr(error, "reason", error)
                error_kind = (
                    "timeout"
                    if isinstance(reason, (socket.timeout, TimeoutError))
                    or "timed out" in str(reason).lower()
                    else "transport_error"
                )
                return self._result(
                    requested_at,
                    method,
                    url,
                    "source_error",
                    None,
                    error_kind,
                    str(reason),
                    None,
                    b"",
                    request_body_sha256,
                )

        raise AssertionError("unreachable HTTP retry state")

    @staticmethod
    def _result(
        requested_at: str,
        method: str,
        url: str,
        status: str,
        http_status: Optional[int],
        error_kind: Optional[str],
        message: Optional[str],
        payload: Any,
        raw: bytes,
        request_body_sha256: Optional[str],
    ) -> dict[str, Any]:
        return {
            "requested_at_utc": requested_at,
            "method": method,
            "url": url,
            "status": status,
            "http_status": http_status,
            "error_kind": error_kind,
            "message": message,
            "payload": payload,
            "response_sha256": hashlib.sha256(raw).hexdigest() if raw else None,
            "request_body_sha256": request_body_sha256,
        }


class FixtureTransport:
    """按 fixture_key 顺序返回固定响应，用于离线测试和可重复验收。"""

    def __init__(
        self, fixtures: dict[str, Any], clock: Callable[[], str] = now_utc
    ):
        self.clock = clock
        self.fixtures: dict[str, list[dict[str, Any]]] = {}
        for key, value in fixtures.items():
            values = value if isinstance(value, list) else [value]
            self.fixtures[key] = [dict(item) for item in values]

    def request_json(
        self,
        fixture_key: str,
        method: str,
        url: str,
        body: Optional[dict[str, Any]] = None,
        body_format: str = "json",
    ) -> dict[str, Any]:
        del body_format
        queue = self.fixtures.get(fixture_key)
        if not queue and ":" in fixture_key:
            queue = self.fixtures.get(fixture_key.split(":", 1)[0])
        if not queue:
            return {
                "requested_at_utc": self.clock(),
                "method": method,
                "url": url,
                "status": "source_error",
                "http_status": None,
                "error_kind": "fixture_missing",
                "message": f"缺少离线响应：{fixture_key}",
                "payload": None,
                "response_sha256": None,
                "request_body_sha256": sha256_json(body) if body else None,
            }
        fixture = queue.pop(0)
        status = fixture.get("status", "success")
        payload = fixture.get("payload")
        raw = canonical_json(payload).encode("utf-8") if payload is not None else b""
        return {
            "requested_at_utc": self.clock(),
            "method": method,
            "url": url,
            "status": status,
            "http_status": fixture.get(
                "http_status",
                200 if status == "success" else 404 if status == "not_found" else None,
            ),
            "error_kind": fixture.get("error_kind"),
            "message": fixture.get("message") or _payload_message(payload),
            "payload": payload,
            "response_sha256": hashlib.sha256(raw).hexdigest() if raw else None,
            "request_body_sha256": sha256_json(body) if body else None,
        }


def _payload_message(payload: Any) -> Optional[str]:
    if not isinstance(payload, dict):
        return None
    for key in ("message", "Message", "response", "Fault"):
        value = payload.get(key)
        if isinstance(value, str):
            return value[:1000]
    fault = payload.get("Fault")
    if isinstance(fault, dict):
        for key in ("Message", "Details", "Code"):
            value = fault.get(key)
            if value:
                return str(value)[:1000]
    return None


def source_log(
    source: str,
    operation: str,
    result: dict[str, Any],
    records_count: int,
    message: Optional[str] = None,
) -> dict[str, Any]:
    return {
        "source": source,
        "operation": operation,
        "requested_at_utc": result["requested_at_utc"],
        "method": result["method"],
        "url": result["url"],
        "status": result["status"],
        "http_status": result.get("http_status"),
        "error_kind": result.get("error_kind"),
        "message": message if message is not None else result.get("message"),
        "records_count": records_count,
        "response_sha256": result.get("response_sha256"),
        "request_body_sha256": result.get("request_body_sha256"),
    }


def fetch_opsin(
    name: str, transport: Any
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    url = f"{OPSIN_BASE}/{urllib.parse.quote(name, safe='')}.json"
    result = transport.request_json("opsin", "GET", url)
    payload = result.get("payload")
    records: list[dict[str, Any]] = []
    message = _payload_message(payload)
    if result["status"] == "success" and isinstance(payload, dict):
        opsin_status = str(payload.get("status") or "").upper()
        if opsin_status in {"SUCCESS", "WARNING"}:
            findings = []
            if opsin_status == "WARNING":
                findings.append(
                    {
                        "code": "R-OPSIN-WARNING",
                        "severity": "review",
                        "message": message or "OPSIN 返回 WARNING。",
                    }
                )
            records.append(
                {
                    "source": "OPSIN",
                    "source_family": "opsin",
                    "source_record_id": None,
                    "match_method": "systematic_name_parser",
                    "title": name,
                    "names": [name],
                    "structure": payload.get("smiles"),
                    "inchi": payload.get("stdinchi"),
                    "inchikey": payload.get("stdinchikey"),
                    "molecular_formula": None,
                    "source_url": url,
                    "raw_record": payload,
                    "record_findings": findings,
                }
            )
        else:
            result = {**result, "status": "not_found", "error_kind": "not_found"}
    elif result["status"] == "not_found" and isinstance(payload, dict):
        message = _payload_message(payload)
    return records, [source_log("OPSIN", "name_to_structure", result, len(records), message)]


def pubchem_request_spec(query: str, input_type: str) -> tuple[str, str, Optional[dict[str, Any]], str]:
    namespace = {
        "name": "name",
        "cas_rn": "name",
        "smiles": "smiles",
        "inchi": "inchi",
        "inchikey": "inchikey",
        "pubchem_cid": "cid",
    }.get(input_type)
    if namespace is None:
        raise ValueError(f"PubChem 不支持 input_type={input_type}")
    if namespace in {"smiles", "inchi"}:
        url = f"{PUBCHEM_BASE}/compound/{namespace}/property/{PUBCHEM_PROPERTIES}/JSON"
        return "POST", url, {namespace: query}, "form"
    encoded = urllib.parse.quote(query, safe="")
    url = (
        f"{PUBCHEM_BASE}/compound/{namespace}/{encoded}/property/"
        f"{PUBCHEM_PROPERTIES}/JSON"
    )
    return "GET", url, None, "json"


def selected_pubchem_record(
    record: dict[str, Any], match_method: str, source_url: str
) -> dict[str, Any]:
    structure = (
        record.get("IsomericSMILES")
        or record.get("SMILES")
        or record.get("CanonicalSMILES")
        or record.get("ConnectivitySMILES")
    )
    return {
        "source": "PubChem",
        "source_family": "pubchem",
        "source_record_id": str(record.get("CID")) if record.get("CID") is not None else None,
        "match_method": match_method,
        "title": record.get("Title"),
        "names": [value for value in (record.get("Title"), record.get("IUPACName")) if value],
        "structure": structure,
        "inchi": record.get("InChI"),
        "inchikey": record.get("InChIKey"),
        "molecular_formula": record.get("MolecularFormula"),
        "source_url": source_url,
        "raw_record": record,
        "record_findings": [],
    }


def fetch_pubchem(
    query: str, input_type: str, transport: Any
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    method, url, body, body_format = pubchem_request_spec(query, input_type)
    result = transport.request_json("pubchem", method, url, body, body_format)
    records: list[dict[str, Any]] = []
    if result["status"] == "success":
        properties = (
            (result.get("payload") or {}).get("PropertyTable", {}).get("Properties", [])
        )
        if isinstance(properties, list):
            records = [
                selected_pubchem_record(item, f"query_{input_type}", url)
                for item in properties
                if isinstance(item, dict)
            ]
        if not records:
            result = {
                **result,
                "status": "not_found",
                "error_kind": "not_found",
                "message": "PubChem 未返回化合物属性记录。",
            }
    return records, [source_log("PubChem", f"lookup_{input_type}", result, len(records))]


def chembl_record_key(record: dict[str, Any]) -> Optional[str]:
    structures = record.get("molecule_structures") or {}
    return structures.get("standard_inchi_key")


def selected_chembl_record(
    record: dict[str, Any], match_method: str, source_url: str, query: Optional[str] = None
) -> dict[str, Any]:
    structures = record.get("molecule_structures") or {}
    properties = record.get("molecule_properties") or {}
    synonyms = record.get("molecule_synonyms") or []
    names = []
    if record.get("pref_name"):
        names.append(record["pref_name"])
    for synonym in synonyms:
        if not isinstance(synonym, dict):
            continue
        value = synonym.get("molecule_synonym")
        if value and (query is None or value.casefold() == query.casefold()):
            names.append(value)
    raw_record = {
        "molecule_chembl_id": record.get("molecule_chembl_id"),
        "pref_name": record.get("pref_name"),
        "molecule_type": record.get("molecule_type"),
        "molecule_hierarchy": record.get("molecule_hierarchy"),
        "molecule_structures": structures,
        "molecule_properties": {
            key: properties.get(key)
            for key in ("full_molformula", "full_mwt")
            if key in properties
        },
        "matched_names": sorted(set(names), key=str.casefold),
    }
    return {
        "source": "ChEMBL",
        "source_family": "chembl",
        "source_record_id": record.get("molecule_chembl_id"),
        "match_method": match_method,
        "title": record.get("pref_name"),
        "names": sorted(set(names), key=str.casefold),
        "structure": structures.get("canonical_smiles"),
        "inchi": structures.get("standard_inchi"),
        "inchikey": structures.get("standard_inchi_key"),
        "molecular_formula": properties.get("full_molformula"),
        "source_url": source_url,
        "raw_record": raw_record,
        "record_findings": [],
    }


def _chembl_list_payload(payload: Any) -> list[dict[str, Any]]:
    if not isinstance(payload, dict):
        return []
    if "molecules" in payload:
        molecules = payload.get("molecules")
        return [item for item in molecules or [] if isinstance(item, dict)]
    if payload.get("molecule_chembl_id"):
        return [payload]
    return []


def fetch_chembl_by_id(
    chembl_id: str, transport: Any
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    normalized = chembl_id.upper()
    url = f"{CHEMBL_BASE}/molecule/{urllib.parse.quote(normalized, safe='')}.json"
    result = transport.request_json("chembl_id", "GET", url)
    raw_records = _chembl_list_payload(result.get("payload"))
    records = [
        selected_chembl_record(record, "chembl_id", url) for record in raw_records
    ]
    if result["status"] == "success" and not records:
        result = {**result, "status": "not_found", "error_kind": "not_found"}
    return records, [source_log("ChEMBL", "lookup_chembl_id", result, len(records))]


def fetch_chembl_by_inchikey(
    inchikey: str, transport: Any
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    query_string = urllib.parse.urlencode(
        {
            "molecule_structures__standard_inchi_key__iexact": inchikey,
            "limit": 20,
        }
    )
    url = f"{CHEMBL_BASE}/molecule.json?{query_string}"
    result = transport.request_json(f"chembl_inchikey:{inchikey}", "GET", url)
    raw_records = _chembl_list_payload(result.get("payload"))
    records = [
        selected_chembl_record(record, "full_inchikey", url)
        for record in raw_records
    ]
    if result["status"] == "success" and not records:
        result = {**result, "status": "not_found", "error_kind": "not_found"}
    return records, [
        source_log("ChEMBL", "lookup_full_inchikey", result, len(records))
    ]


def fetch_chembl_by_name(
    name: str, transport: Any
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    records_by_id: dict[str, dict[str, Any]] = {}
    logs: list[dict[str, Any]] = []
    operations = [
        (
            "chembl_pref_name",
            "pref_name__iexact",
            "exact_preferred_name",
        ),
        (
            "chembl_synonym",
            "molecule_synonyms__molecule_synonym__iexact",
            "exact_synonym",
        ),
    ]
    for fixture_key, field, match_method in operations:
        query_string = urllib.parse.urlencode({field: name, "limit": 50})
        url = f"{CHEMBL_BASE}/molecule.json?{query_string}"
        result = transport.request_json(fixture_key, "GET", url)
        raw_records = _chembl_list_payload(result.get("payload"))
        selected = [
            selected_chembl_record(record, match_method, url, name)
            for record in raw_records
        ]
        for record in selected:
            identifier = record.get("source_record_id") or sha256_json(record)
            existing = records_by_id.get(identifier)
            if existing is None:
                records_by_id[identifier] = record
            else:
                existing["names"] = sorted(
                    set(existing["names"] + record["names"]), key=str.casefold
                )
                existing["match_method"] = "exact_preferred_name_and_synonym"
        if result["status"] == "success" and not selected:
            result = {**result, "status": "not_found", "error_kind": "not_found"}
        logs.append(
            source_log("ChEMBL", f"lookup_{match_method}", result, len(selected))
        )
    return list(records_by_id.values()), logs


def selected_unichem_record(
    compound: dict[str, Any], source_url: str
) -> dict[str, Any]:
    inchi_data = compound.get("inchi") or {}
    sources = compound.get("sources") or []
    selected_sources = [
        {
            "source": source.get("shortName"),
            "source_record_id": source.get("compoundId"),
            "source_url": source.get("url"),
        }
        for source in sources
        if isinstance(source, dict)
        and source.get("shortName") in {"pubchem", "chembl"}
    ]
    raw_record = {
        "uci": compound.get("uci"),
        "standardInchiKey": compound.get("standardInchiKey"),
        "inchi": inchi_data,
        "selected_source_mappings": selected_sources,
        "source_mapping_count": len(sources),
    }
    return {
        "source": "UniChem",
        "source_family": "unichem",
        "source_record_id": str(compound.get("uci")) if compound.get("uci") is not None else None,
        "match_method": "exact_full_inchikey",
        "title": None,
        "names": [],
        "structure": None,
        "inchi": inchi_data.get("inchi"),
        "inchikey": compound.get("standardInchiKey"),
        "molecular_formula": inchi_data.get("formula"),
        "source_url": source_url,
        "raw_record": raw_record,
        "record_findings": [],
    }


def fetch_unichem_exact(
    inchikey: str, transport: Any
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    url = f"{UNICHEM_BASE}/compounds"
    body = {"compound": inchikey, "type": "inchikey"}
    result = transport.request_json(
        f"unichem_exact:{inchikey}", "POST", url, body, "json"
    )
    payload = result.get("payload")
    records = []
    if result["status"] == "success" and isinstance(payload, dict):
        compounds = payload.get("compounds") or []
        records = [
            selected_unichem_record(compound, url)
            for compound in compounds
            if isinstance(compound, dict)
        ]
        if not records:
            result = {
                **result,
                "status": "not_found",
                "error_kind": "not_found",
                "message": "UniChem exact 未返回匹配。",
            }
    return records, [source_log("UniChem", "exact_full_inchikey", result, len(records))]


def fetch_unichem_connectivity(
    inchikey: str, transport: Any
) -> tuple[Optional[dict[str, Any]], list[dict[str, Any]]]:
    url = f"{UNICHEM_BASE}/connectivity"
    body = {
        "compound": inchikey,
        "type": "inchikey",
        "searchComponents": True,
    }
    result = transport.request_json(
        f"unichem_connectivity:{inchikey}", "POST", url, body, "json"
    )
    payload = result.get("payload")
    summary = None
    if result["status"] == "success" and isinstance(payload, dict):
        sources = [item for item in payload.get("sources") or [] if isinstance(item, dict)]
        type_counts: dict[str, int] = {}
        layer_difference_counts: dict[str, int] = {}
        for source in sources:
            search_type = str(source.get("typeOfSearch") or "unknown")
            type_counts[search_type] = type_counts.get(search_type, 0) + 1
            comparison = source.get("comparison") or {}
            for layer, same in comparison.items():
                if same is False:
                    layer_difference_counts[layer] = (
                        layer_difference_counts.get(layer, 0) + 1
                    )
        selected_sources = [
            {
                "source": source.get("shortName"),
                "source_record_id": source.get("compoundId"),
                "source_url": source.get("url"),
                "type_of_search": source.get("typeOfSearch"),
                "searched_component_index": source.get("searchedInchiPos"),
                "layer_comparison": source.get("comparison"),
            }
            for source in sources
            if source.get("shortName") in {"pubchem", "chembl"}
        ][:100]
        summary = {
            "query_inchikey": inchikey,
            "searched_compound": payload.get("searchedCompound"),
            "total_compounds": payload.get("totalCompounds"),
            "total_sources": payload.get("totalSources"),
            "search_type_counts": type_counts,
            "layer_difference_counts": layer_difference_counts,
            "selected_source_records": selected_sources,
            "selected_source_records_truncated": len(selected_sources) >= 100,
        }
    return summary, [
        source_log(
            "UniChem",
            "connectivity_related_forms",
            result,
            len((summary or {}).get("selected_source_records", [])),
        )
    ]


def fetch_source_metadata(
    enabled_sources: set[str], transport: Any
) -> dict[str, Any]:
    metadata: dict[str, Any] = {
        "OPSIN": {
            "enabled": "opsin" in enabled_sources,
            "service_version": "not_exposed_by_api",
            "reviewed_release": "2.9.0",
            "documentation_url": "https://github.com/dan2097/opsin/releases/tag/2.9.0",
        },
        "PubChem": {
            "enabled": "pubchem" in enabled_sources,
            "service_version": "not_exposed_by_api",
            "api": "PUG REST",
            "documentation_url": "https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest",
        },
        "ChEMBL": {
            "enabled": "chembl" in enabled_sources,
            "database_version": None,
            "release_date": None,
            "api_status": None,
            "documentation_url": (
                "https://chembl.gitbook.io/chembl-interface-documentation/"
                "web-services/chembl-data-web-services"
            ),
        },
        "UniChem": {
            "enabled": "unichem" in enabled_sources,
            "service_version": "not_exposed_by_api",
            "api": "POST API v1 (2.0 documentation)",
            "documentation_url": "https://chembl.gitbook.io/unichem/api",
        },
    }
    if "chembl" not in enabled_sources:
        return metadata

    url = f"{CHEMBL_BASE}/status.json"
    result = transport.request_json("chembl_status", "GET", url)
    payload = result.get("payload")
    if result["status"] == "success" and isinstance(payload, dict):
        metadata["ChEMBL"].update(
            {
                "database_version": payload.get("chembl_db_version"),
                "release_date": payload.get("chembl_release_date"),
                "api_status": payload.get("status"),
                "database_counts": {
                    key: payload.get(key)
                    for key in (
                        "activities",
                        "compound_records",
                        "disinct_compounds",
                        "publications",
                        "targets",
                    )
                    if key in payload
                },
            }
        )
    metadata["ChEMBL"]["status_query"] = source_log(
        "ChEMBL", "runtime_status", result, 1 if isinstance(payload, dict) else 0
    )
    return metadata


def normalize_source_record(
    record: dict[str, Any], toolkit: dict[str, Any]
) -> tuple[Optional[dict[str, Any]], Optional[dict[str, Any]]]:
    structure = record.get("structure")
    source_inchi = record.get("inchi")
    source_key = record.get("inchikey")
    molecule = None
    parse_error = None
    if source_inchi:
        molecule, parse_error = parse_structure(source_inchi, "inchi", toolkit)
    if molecule is None and structure:
        molecule, parse_error = parse_structure(structure, "smiles", toolkit)

    derived = None
    findings = list(record.get("record_findings") or [])
    if molecule is not None:
        derived = structure_identifiers(molecule, toolkit)
    if derived and source_key and derived["inchikey"] != source_key:
        conflict = {
            "code": "E-SOURCE-INCHIKEY-MISMATCH",
            "severity": "error",
            "message": (
                f"{record['source']} 提供的 InChIKey 与本地从其结构重算的完整 "
                "InChIKey 不一致。"
            ),
            "source": record["source"],
            "source_record_id": record.get("source_record_id"),
            "source_inchikey": source_key,
            "derived_inchikey": derived["inchikey"],
        }
        findings.append(conflict)
        return None, {**record, "record_findings": findings, "parse_error": parse_error}

    final_key = (derived or {}).get("inchikey") or source_key
    if not final_key:
        unresolved = {
            **record,
            "record_findings": findings
            + [
                {
                    "code": "R-SOURCE-RECORD-NO-COMPLETE-STRUCTURE",
                    "severity": "review",
                    "message": "来源记录没有可核对的完整结构或 InChIKey。",
                }
            ],
            "parse_error": parse_error,
        }
        return None, unresolved

    if derived is None:
        findings.append(
            {
                "code": "R-KEY-WITHOUT-LOCAL-STRUCTURE",
                "severity": "review",
                "message": "来源只提供完整 InChIKey，未能在本地重建结构。",
            }
        )
    normalized = {
        **record,
        "canonical_smiles": (derived or {}).get("canonical_smiles"),
        "normalized_inchi": (derived or {}).get("inchi") or source_inchi,
        "normalized_inchikey": final_key,
        "connectivity_block": final_key[:14] if INCHIKEY_RE.fullmatch(final_key) else None,
        "normalized_formula": (derived or {}).get("molecular_formula")
        or record.get("molecular_formula"),
        "component_count": (derived or {}).get("component_count"),
        "unassigned_stereo": (derived or {}).get("unassigned_stereo") or [],
        "record_findings": findings,
    }
    if normalized["component_count"] and normalized["component_count"] > 1:
        normalized["record_findings"].append(
            {
                "code": "R-MULTICOMPONENT-CANDIDATE",
                "severity": "review",
                "message": "候选结构包含多个组分，可能是盐、溶剂化物、配合物或混合物。",
            }
        )
    if normalized["unassigned_stereo"]:
        normalized["record_findings"].append(
            {
                "code": "R-UNSPECIFIED-STEREO",
                "severity": "review",
                "message": "候选结构存在未指定立体中心。",
            }
        )
    return normalized, None


def aggregate_candidates(
    records: Sequence[dict[str, Any]], toolkit: dict[str, Any]
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    groups: dict[str, dict[str, Any]] = {}
    unresolved: list[dict[str, Any]] = []
    integrity_conflicts: list[dict[str, Any]] = []
    for record in records:
        normalized, problem = normalize_source_record(record, toolkit)
        if normalized is None:
            if problem and any(
                item.get("code") == "E-SOURCE-INCHIKEY-MISMATCH"
                for item in problem.get("record_findings", [])
            ):
                integrity_conflicts.append(problem)
            elif problem:
                unresolved.append(problem)
            continue
        key = normalized["normalized_inchikey"]
        group = groups.setdefault(
            key,
            {
                "canonical_smiles": normalized.get("canonical_smiles"),
                "inchi": normalized.get("normalized_inchi"),
                "inchikey": key,
                "connectivity_block": normalized.get("connectivity_block"),
                "molecular_formula": normalized.get("normalized_formula"),
                "component_count": normalized.get("component_count"),
                "names": set(),
                "evidence": [],
                "quality_findings": [],
                "comparison_view": None,
            },
        )
        if group["canonical_smiles"] is None and normalized.get("canonical_smiles"):
            group["canonical_smiles"] = normalized["canonical_smiles"]
        if group["inchi"] is None and normalized.get("normalized_inchi"):
            group["inchi"] = normalized["normalized_inchi"]
        if group["molecular_formula"] is None and normalized.get("normalized_formula"):
            group["molecular_formula"] = normalized["normalized_formula"]
        group["names"].update(normalized.get("names") or [])
        if normalized.get("title"):
            group["names"].add(normalized["title"])
        evidence = {
            "source": normalized["source"],
            "source_family": normalized["source_family"],
            "source_record_id": normalized.get("source_record_id"),
            "source_url": normalized.get("source_url"),
            "match_method": normalized.get("match_method"),
            "source_structure": normalized.get("structure"),
            "source_inchi": normalized.get("inchi"),
            "source_inchikey": normalized.get("inchikey"),
            "raw_record": normalized.get("raw_record"),
        }
        evidence_key = (
            evidence["source_family"],
            evidence["source_record_id"],
            evidence["match_method"],
            evidence["source_inchikey"],
        )
        existing_keys = {
            (
                item["source_family"],
                item["source_record_id"],
                item["match_method"],
                item["source_inchikey"],
            )
            for item in group["evidence"]
        }
        if evidence_key not in existing_keys:
            group["evidence"].append(evidence)
        for finding in normalized.get("record_findings") or []:
            finding_key = (finding.get("code"), finding.get("message"))
            existing_findings = {
                (item.get("code"), item.get("message"))
                for item in group["quality_findings"]
            }
            if finding_key not in existing_findings:
                group["quality_findings"].append(finding)

    sorted_groups = sorted(
        groups.values(),
        key=lambda item: (
            item["inchikey"] or "",
            item["canonical_smiles"] or "",
        ),
    )
    candidates = []
    for index, group in enumerate(sorted_groups, 1):
        group["candidate_id"] = f"candidate-{index:03d}"
        group["names"] = sorted(group["names"], key=str.casefold)
        group["evidence"] = sorted(
            group["evidence"],
            key=lambda item: (
                item["source_family"],
                item.get("source_record_id") or "",
                item.get("match_method") or "",
            ),
        )
        group["source_families"] = sorted(
            {item["source_family"] for item in group["evidence"]}
        )
        candidates.append(group)
    return candidates, unresolved, integrity_conflicts


def default_standardizer_path() -> Path:
    return (
        Path(__file__).resolve().parents[2]
        / "standardize-chemical-structures"
        / "scripts"
        / "standardize_structures.py"
    )


def apply_standardization_views(
    candidates: list[dict[str, Any]],
    standardizer_script: Optional[Path],
    profile: str,
    generated_at_utc: str,
) -> dict[str, Any]:
    eligible = [candidate for candidate in candidates if candidate["canonical_smiles"]]
    if not eligible:
        return {
            "status": "not_run",
            "reason": "没有可交给结构标准化 Skill 的候选结构。",
            "target_skill": "standardize-chemical-structures",
            "profile": profile,
        }
    if standardizer_script is None:
        return {
            "status": "not_run",
            "reason": "调用方显式禁用了结构标准化交接。",
            "target_skill": "standardize-chemical-structures",
            "profile": profile,
        }
    if not standardizer_script.exists():
        return {
            "status": "not_run",
            "reason": f"未找到结构标准化脚本：{standardizer_script}",
            "target_skill": "standardize-chemical-structures",
            "profile": profile,
        }

    spec = importlib.util.spec_from_file_location(
        "_identity_standardizer", standardizer_script
    )
    if spec is None or spec.loader is None:
        return {
            "status": "error",
            "reason": "无法加载结构标准化脚本。",
            "target_skill": "standardize-chemical-structures",
            "profile": profile,
        }
    module = importlib.util.module_from_spec(spec)
    try:
        spec.loader.exec_module(module)
        input_records = [
            {
                "id": candidate["candidate_id"],
                "record_index": index,
                "source": "resolve-chemical-identities",
                "input_format": "smiles",
                "original_structure": candidate["canonical_smiles"],
            }
            for index, candidate in enumerate(eligible, 1)
        ]
        document = module.process_records(
            input_records,
            profile,
            provenance=[
                {
                    "source": "resolve-chemical-identities",
                    "purpose": "derived_comparison_view_only",
                }
            ],
            generated_at_utc=generated_at_utc,
        )
    except Exception as error:
        return {
            "status": "error",
            "reason": str(error),
            "target_skill": "standardize-chemical-structures",
            "profile": profile,
        }

    by_id = {record["id"]: record for record in document["records"]}
    for candidate in eligible:
        record = by_id[candidate["candidate_id"]]
        candidate["comparison_view"] = {
            "status": record["standardization_status"],
            "profile": profile,
            "standardized_structure": record["standardized_structure"],
            "parent_structure": record["parent_structure"],
            "standardized_inchikey": record["inchikey"],
            "parent_inchikey": record["parent_inchikey"],
            "disposition": record["disposition"],
            "finding_codes": [
                finding["code"] for finding in record["qc_findings"]
            ],
        }
    return {
        "status": "completed",
        "target_skill": "standardize-chemical-structures",
        "profile": profile,
        "tool_versions": document["tool_versions"],
        "notice": (
            "comparison_view 是派生比较口径，不覆盖来源结构；"
            "同一 parent 不代表同一物理样品。"
        ),
    }


def aggregate_retrieval_status(source_queries: Sequence[dict[str, Any]]) -> str:
    if not source_queries:
        return "not_run"
    statuses = [item["status"] for item in source_queries]
    successes = statuses.count("success")
    errors = statuses.count("source_error")
    if successes and errors:
        return "partial"
    if errors and not successes:
        return "source_error"
    if successes:
        return "completed"
    if all(status == "not_found" for status in statuses):
        return "not_found"
    return "source_error"


def has_review_findings(candidate: dict[str, Any]) -> bool:
    if any(
        finding.get("severity") in {"review", "error"}
        for finding in candidate.get("quality_findings", [])
    ):
        return True
    comparison = candidate.get("comparison_view") or {}
    return comparison.get("disposition") in {"review_required", "rejected"}


def determine_alignment(
    validated: dict[str, Any],
    candidates: Sequence[dict[str, Any]],
    unresolved: Sequence[dict[str, Any]],
    integrity_conflicts: Sequence[dict[str, Any]],
    retrieval_status: str,
) -> tuple[str, str, list[dict[str, Any]]]:
    questions: list[dict[str, Any]] = []
    input_type = validated["detected_input_type"]

    if validated["input_status"] == "invalid_input":
        return "not_assessed", "rejected", questions
    if integrity_conflicts:
        questions.append(
            {
                "code": "Q-SOURCE-INTEGRITY-CONFLICT",
                "question": "来源结构与来源 InChIKey 冲突，应以哪个经人工核验的记录为准？",
            }
        )
        return "conflict", "review_required", questions
    if not candidates:
        if retrieval_status == "not_found":
            return "not_assessed", "rejected", questions
        questions.append(
            {
                "code": "Q-NO-RESOLVED-STRUCTURE",
                "question": "当前没有可核对的完整结构；是否能提供结构、稳定 ID 或更多上下文？",
            }
        )
        return "not_assessed", "review_required", questions

    if len(candidates) > 1:
        parent_keys = [
            (candidate.get("comparison_view") or {}).get("parent_inchikey")
            for candidate in candidates
        ]
        same_parent = all(parent_keys) and len(set(parent_keys)) == 1
        if same_parent:
            questions.append(
                {
                    "code": "Q-RELATED-FORM-SELECTION",
                    "question": (
                        "候选共享派生 parent，但完整结构不同；"
                        "请确认需要哪一种盐型、质子化形式、立体或同位素形式。"
                    ),
                }
            )
            return "related_forms", "review_required", questions
        if input_type in {"name"}:
            questions.append(
                {
                    "code": "Q-AMBIGUOUS-NAME",
                    "question": "该名称对应多个完整结构；请提供盐型、立体、用途或稳定标识符。",
                }
            )
            return "ambiguous", "review_required", questions
        questions.append(
            {
                "code": "Q-STABLE-ID-CONFLICT",
                "question": "同一结构或稳定标识符返回多个不一致的完整结构，请人工核对来源。",
            }
        )
        return "conflict", "review_required", questions

    candidate = candidates[0]
    if unresolved and input_type == "name":
        questions.append(
            {
                "code": "Q-STRUCTURELESS-NAME-RECORD",
                "question": (
                    "至少一个精确名称记录没有完整结构，可能代表家族或未指定形式；"
                    "请确认具体化学形式。"
                ),
            }
        )
        return "ambiguous", "review_required", questions

    evidence_families = set(candidate["source_families"])
    external_families = evidence_families - {"local_input", "unichem"}
    if input_type == "name":
        if len(external_families) < 2:
            questions.append(
                {
                    "code": "Q-SINGLE-SOURCE-NAME",
                    "question": "名称目前只有一个独立证据源支持，是否能提供更多上下文或稳定 ID？",
                }
            )
            return "ambiguous", "review_required", questions
        alignment = "exact"
    elif input_type in {"smiles", "inchi"} and not external_families:
        alignment = "not_assessed"
    else:
        alignment = "exact"

    if retrieval_status in {"partial", "source_error"}:
        questions.append(
            {
                "code": "Q-RETRY-SOURCE-ERROR",
                "question": "至少一个来源查询失败；是否在服务恢复后重试再确认记录对齐？",
            }
        )
        return alignment, "review_required", questions
    if has_review_findings(candidate):
        questions.append(
            {
                "code": "Q-CANDIDATE-QUALITY-REVIEW",
                "question": "候选包含多组分、未指定立体或标准化复核项，请确认是否适合目标用途。",
            }
        )
        return alignment, "review_required", questions
    if alignment == "not_assessed" and input_type in {"smiles", "inchi"}:
        return alignment, "ready_for_standardization", questions
    if alignment == "exact":
        return alignment, "ready_for_standardization", questions
    return alignment, "review_required", questions


def build_handoff(
    validated: dict[str, Any],
    candidates: Sequence[dict[str, Any]],
    alignment: str,
    disposition: str,
) -> dict[str, Any]:
    if disposition != "ready_for_standardization" or len(candidates) != 1:
        return {
            "status": "blocked_pending_resolution",
            "target_skill": "standardize-chemical-structures",
            "records": [],
            "reason": "只有唯一候选且当前处置为 ready_for_standardization 时才生成下游输入。",
        }
    candidate = candidates[0]
    structure = candidate.get("canonical_smiles")
    if not structure:
        return {
            "status": "blocked_missing_structure",
            "target_skill": "standardize-chemical-structures",
            "records": [],
            "reason": "唯一候选没有可交付的结构。",
        }
    return {
        "status": "ready",
        "target_skill": "standardize-chemical-structures",
        "records": [
            {
                "id": validated["id"],
                "structure": structure,
                "source_candidate_id": candidate["candidate_id"],
                "source_inchikey": candidate["inchikey"],
            }
        ],
        "alignment_scope": (
            "input_structure_only"
            if alignment == "not_assessed"
            else "database_records_only"
        ),
        "notice": (
            "该交接不确认物理样品身份；第一个 Skill 仍须保留此处来源结构和候选证据。"
        ),
    }


def resolve_one(
    item: dict[str, Any],
    toolkit: dict[str, Any],
    transport: Any,
    enabled_sources: set[str],
    include_related: bool,
    standardizer_script: Optional[Path],
    standardization_profile: str,
    generated_at_utc: str,
) -> dict[str, Any]:
    validated = validate_request(item, toolkit)
    source_records: list[dict[str, Any]] = []
    source_queries: list[dict[str, Any]] = []
    relationship_evidence: list[dict[str, Any]] = []

    if validated.get("local_record"):
        source_records.append(validated["local_record"])
    if validated["input_status"] == "invalid_input":
        return {
            "request": {
                key: validated.get(key)
                for key in (
                    "id",
                    "query",
                    "normalized_query",
                    "requested_input_type",
                    "detected_input_type",
                    "context",
                    "expected_form",
                )
            },
            "input_status": "invalid_input",
            "retrieval_status": "not_run",
            "record_alignment_status": "not_assessed",
            "record_alignment_scope": "database_records_only",
            "sample_identity_status": "not_assessed",
            "disposition": "rejected",
            "candidates": [],
            "unresolved_source_records": [],
            "source_record_conflicts": [],
            "source_queries": [],
            "relationship_evidence": [],
            "confirmation_questions": [],
            "findings": validated["findings"],
            "standardization_comparison": {
                "status": "not_run",
                "reason": "输入无效。",
                "target_skill": "standardize-chemical-structures",
            },
            "standardization_handoff": {
                "status": "blocked_invalid_input",
                "target_skill": "standardize-chemical-structures",
                "records": [],
            },
        }

    query = validated["normalized_query"]
    input_type = validated["detected_input_type"]

    if "opsin" in enabled_sources and input_type == "name":
        records, logs = fetch_opsin(query, transport)
        source_records.extend(records)
        source_queries.extend(logs)

    if "pubchem" in enabled_sources and input_type != "chembl_id":
        records, logs = fetch_pubchem(query, input_type, transport)
        source_records.extend(records)
        source_queries.extend(logs)

    if "chembl" in enabled_sources:
        if input_type == "chembl_id":
            records, logs = fetch_chembl_by_id(query, transport)
            source_records.extend(records)
            source_queries.extend(logs)
        elif input_type in {"name", "cas_rn"}:
            records, logs = fetch_chembl_by_name(query, transport)
            source_records.extend(records)
            source_queries.extend(logs)
        else:
            preliminary, _, _ = aggregate_candidates(source_records, toolkit)
            keys = sorted(
                {candidate["inchikey"] for candidate in preliminary if candidate["inchikey"]}
            )
            for key in keys:
                records, logs = fetch_chembl_by_inchikey(key, transport)
                source_records.extend(records)
                source_queries.extend(logs)

    if input_type == "chembl_id" and "pubchem" in enabled_sources:
        preliminary, _, _ = aggregate_candidates(source_records, toolkit)
        keys = sorted(
            {candidate["inchikey"] for candidate in preliminary if candidate["inchikey"]}
        )
        for key in keys:
            records, logs = fetch_pubchem(key, "inchikey", transport)
            source_records.extend(records)
            source_queries.extend(logs)

    if "unichem" in enabled_sources:
        preliminary, _, _ = aggregate_candidates(source_records, toolkit)
        keys = sorted(
            {candidate["inchikey"] for candidate in preliminary if candidate["inchikey"]}
        )
        for key in keys:
            records, logs = fetch_unichem_exact(key, transport)
            source_records.extend(records)
            source_queries.extend(logs)
            if include_related:
                summary, relation_logs = fetch_unichem_connectivity(key, transport)
                source_queries.extend(relation_logs)
                if summary:
                    relationship_evidence.append(summary)

    candidates, unresolved, integrity_conflicts = aggregate_candidates(
        source_records, toolkit
    )
    standardization_comparison = apply_standardization_views(
        candidates,
        standardizer_script,
        standardization_profile,
        generated_at_utc,
    )
    retrieval_status = aggregate_retrieval_status(source_queries)
    alignment, disposition, questions = determine_alignment(
        validated,
        candidates,
        unresolved,
        integrity_conflicts,
        retrieval_status,
    )

    if validated.get("expected_form") and disposition == "ready_for_standardization":
        # expected_form is user context, not proof. A mismatch must be reviewed manually.
        expected_form = str(validated["expected_form"]).strip().lower()
        candidate = candidates[0]
        if expected_form in {"neutral", "single_component"} and (
            candidate.get("component_count") or 1
        ) > 1:
            disposition = "review_required"
            questions.append(
                {
                    "code": "Q-EXPECTED-FORM-MISMATCH",
                    "question": "用户期望单一中性形式，但候选包含多个组分；请确认具体形式。",
                }
            )

    handoff = build_handoff(validated, candidates, alignment, disposition)
    findings = list(validated["findings"])
    if retrieval_status == "partial":
        findings.append(
            {
                "code": "R-PARTIAL-SOURCE-RETRIEVAL",
                "severity": "review",
                "message": "至少一个来源成功且至少一个来源故障；不能把故障当成无记录。",
            }
        )
    elif retrieval_status == "source_error":
        findings.append(
            {
                "code": "R-SOURCE-ERROR",
                "severity": "review",
                "message": "来源查询失败；当前结果不能表述为 not_found。",
            }
        )
    elif retrieval_status == "not_found":
        findings.append(
            {
                "code": "E-NOT-FOUND",
                "severity": "error",
                "message": "已查询来源均明确未返回记录。",
            }
        )

    return {
        "request": {
            key: validated.get(key)
            for key in (
                "id",
                "query",
                "normalized_query",
                "requested_input_type",
                "detected_input_type",
                "context",
                "expected_form",
            )
        },
        "input_status": validated["input_status"],
        "retrieval_status": retrieval_status,
        "record_alignment_status": alignment,
        "record_alignment_scope": "database_records_only",
        "sample_identity_status": "not_assessed",
        "disposition": disposition,
        "candidates": candidates,
        "unresolved_source_records": unresolved,
        "source_record_conflicts": integrity_conflicts,
        "source_queries": source_queries,
        "relationship_evidence": relationship_evidence,
        "confirmation_questions": questions,
        "findings": findings,
        "standardization_comparison": standardization_comparison,
        "standardization_handoff": handoff,
    }


def build_cross_query_relationships(
    resolutions: Sequence[dict[str, Any]]
) -> list[dict[str, Any]]:
    relationships = []
    for left_index, left in enumerate(resolutions):
        if len(left["candidates"]) != 1:
            continue
        left_candidate = left["candidates"][0]
        for right in resolutions[left_index + 1 :]:
            if len(right["candidates"]) != 1:
                continue
            right_candidate = right["candidates"][0]
            left_key = left_candidate.get("inchikey")
            right_key = right_candidate.get("inchikey")
            left_parent = (left_candidate.get("comparison_view") or {}).get(
                "parent_inchikey"
            )
            right_parent = (right_candidate.get("comparison_view") or {}).get(
                "parent_inchikey"
            )
            if left_key and left_key == right_key:
                relationship = "exact"
                explanation = "两个查询解析到相同完整 InChIKey，仅表示数字记录结构一致。"
            elif left_parent and left_parent == right_parent:
                relationship = "related_forms"
                explanation = (
                    "两个完整结构不同但共享派生 parent；"
                    "可能是盐型或其他相关形式，不是同一物理样品。"
                )
            else:
                relationship = "different_or_unresolved"
                explanation = "当前确定性规则未证明两个查询为相同记录或相关形式。"
            relationships.append(
                {
                    "left_request_id": left["request"]["id"],
                    "right_request_id": right["request"]["id"],
                    "relationship": relationship,
                    "left_inchikey": left_key,
                    "right_inchikey": right_key,
                    "left_parent_inchikey": left_parent,
                    "right_parent_inchikey": right_parent,
                    "explanation": explanation,
                }
            )
    return relationships


def process_requests(
    requests: Sequence[dict[str, Any]],
    *,
    transport: Any,
    enabled_sources: Iterable[str] = DEFAULT_SOURCES,
    include_related: bool = False,
    use_standardizer: bool = True,
    standardizer_script: Optional[Path] = None,
    standardization_profile: str = "chembl-pipeline",
    generated_at_utc: Optional[str] = None,
) -> dict[str, Any]:
    if not requests:
        raise InputFailure("至少需要一个 query。")
    toolkit = load_toolkit()
    source_set = {source.strip().lower() for source in enabled_sources if source}
    unsupported = source_set - SUPPORTED_SOURCES
    if unsupported:
        raise InputFailure(f"不支持的来源：{', '.join(sorted(unsupported))}")
    generated_at = generated_at_utc or now_utc()
    if use_standardizer and standardizer_script is None:
        standardizer_script = default_standardizer_path()
    elif not use_standardizer:
        standardizer_script = None
    source_metadata = fetch_source_metadata(source_set, transport)
    resolutions = [
        resolve_one(
            dict(item),
            toolkit,
            transport,
            source_set,
            include_related,
            standardizer_script,
            standardization_profile,
            generated_at,
        )
        for item in requests
    ]
    counts = {
        status: sum(
            resolution["disposition"] == status for resolution in resolutions
        )
        for status in ("ready_for_standardization", "review_required", "rejected")
    }
    document = {
        "schema_version": SCHEMA_VERSION,
        "workflow": WORKFLOW,
        "generated_at_utc": generated_at,
        "tool_versions": toolkit_versions(toolkit),
        "source_metadata": source_metadata,
        "options": {
            "enabled_sources": sorted(source_set),
            "include_unichem_connectivity": include_related,
            "standardization_profile": standardization_profile,
            "use_standardizer": use_standardizer,
            "standardizer_script": (
                str(standardizer_script) if standardizer_script else None
            ),
            "no_model_generated_structures": True,
            "automatic_tie_breaking": False,
        },
        "input_summary": {
            "total_requests": len(resolutions),
            **counts,
        },
        "resolutions": resolutions,
        "cross_query_relationships": build_cross_query_relationships(resolutions),
        "notices": [
            "record_alignment_status 只描述输入与数字来源记录的结构关系，不确认用户实物样品。",
            "sample_identity_status 默认 not_assessed；数据库数量、评分或首选记录不能自动升级该状态。",
            "comparison_view 和 parent 是派生比较口径，不覆盖来源结构，也不表示同一物理样品。",
            "本工作流不判断活性、药效、毒性、可合成性或实验安全。",
        ],
    }
    document["result_fingerprint"] = output_fingerprint(document)
    serialized = json.dumps(document, ensure_ascii=False)
    if SECRET_RE.search(serialized):
        raise RuntimeError("输出中检测到疑似凭证，已停止写出。")
    return document


def load_request_file(path: Path) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise InputFailure(f"无法读取请求文件：{error}") from error
    if not isinstance(payload, dict):
        raise InputFailure("请求文件顶层必须是 JSON object。")
    requests = payload.get("requests")
    if requests is None and "query" in payload:
        requests = [payload]
    if not isinstance(requests, list) or not all(
        isinstance(item, dict) for item in requests
    ):
        raise InputFailure("请求文件必须包含 requests 数组或单个 query。")
    options = payload.get("options") or {}
    if not isinstance(options, dict):
        raise InputFailure("options 必须是 JSON object。")
    return requests, options


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--query", action="append", default=[])
    parser.add_argument("--request", type=Path, help="单条或批量 JSON 请求文件")
    parser.add_argument("--input-type", default="auto", choices=sorted(INPUT_TYPES))
    parser.add_argument("--context")
    parser.add_argument("--expected-form")
    parser.add_argument(
        "--sources",
        default=",".join(DEFAULT_SOURCES),
        help="逗号分隔：opsin,pubchem,chembl,unichem；空字符串表示不联网",
    )
    parser.add_argument(
        "--include-related",
        action="store_true",
        help="额外调用 UniChem connectivity，记录相关形式证据",
    )
    parser.add_argument(
        "--standardizer-script",
        type=Path,
        help="第一个 Skill 的 standardize_structures.py；默认自动寻找同级 Skill",
    )
    parser.add_argument(
        "--no-standardizer",
        action="store_true",
        help="不生成第一个 Skill 的派生 comparison_view",
    )
    parser.add_argument(
        "--standardization-profile",
        default="chembl-pipeline",
        choices=["rdkit-basic", "chembl-pipeline"],
    )
    parser.add_argument("--timeout", type=int, default=DEFAULT_TIMEOUT)
    parser.add_argument("--retries", type=int, default=DEFAULT_RETRIES)
    parser.add_argument(
        "--fixture-responses",
        type=Path,
        help="仅用于离线测试的固定来源响应 JSON",
    )
    parser.add_argument("--generated-at", help="固定 UTC 时间，仅用于可重复验收")
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        file_options: dict[str, Any] = {}
        if args.request:
            requests, file_options = load_request_file(args.request)
            if args.query:
                raise InputFailure("--request 与 --query 不能同时使用。")
        else:
            requests = [
                {
                    "id": f"query-{index}",
                    "query": query,
                    "input_type": args.input_type,
                    "context": args.context,
                    "expected_form": args.expected_form,
                }
                for index, query in enumerate(args.query, 1)
            ]
        if not requests:
            raise InputFailure("请提供 --query 或 --request。")

        generated_at = args.generated_at or now_utc()

        def clock() -> str:
            return generated_at

        if args.fixture_responses:
            fixtures = json.loads(args.fixture_responses.read_text(encoding="utf-8"))
            if not isinstance(fixtures, dict):
                raise InputFailure("fixture-responses 顶层必须是 JSON object。")
            transport = FixtureTransport(fixtures, clock=clock)
        else:
            transport = HttpTransport(
                timeout=args.timeout,
                retries=args.retries,
                clock=clock,
            )

        sources_value = file_options.get("sources", args.sources)
        if isinstance(sources_value, str):
            sources = [
                source.strip()
                for source in sources_value.split(",")
                if source.strip()
            ]
        elif isinstance(sources_value, list):
            sources = [str(source) for source in sources_value]
        else:
            raise InputFailure("sources 必须是逗号分隔字符串或数组。")
        include_related = bool(
            file_options.get("include_related", args.include_related)
        )
        standardizer_script = args.standardizer_script or default_standardizer_path()
        document = process_requests(
            requests,
            transport=transport,
            enabled_sources=sources,
            include_related=include_related,
            use_standardizer=not args.no_standardizer,
            standardizer_script=standardizer_script,
            standardization_profile=file_options.get(
                "standardization_profile", args.standardization_profile
            ),
            generated_at_utc=generated_at,
        )
    except (
        DependencyFailure,
        InputFailure,
        OSError,
        ValueError,
        json.JSONDecodeError,
    ) as error:
        sys.stderr.write(f"error: {error}\n")
        return 3

    serialized = json.dumps(document, ensure_ascii=False, indent=2) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(serialized, encoding="utf-8")
    else:
        sys.stdout.write(serialized)
    return 2 if document["input_summary"]["rejected"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
