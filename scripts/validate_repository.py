#!/usr/bin/env python3
from __future__ import annotations

import re
import sys
from pathlib import Path
from urllib.parse import unquote, urlsplit

import yaml


ROOT = Path(__file__).resolve().parents[1]
SKILLS_ROOT = ROOT / "skills"

SKILLS = (
    "resolve-chemical-identities",
    "standardize-chemical-structures",
    "compute-molecular-features",
    "search-and-curate-chemical-libraries",
    "curate-reactions",
    "search-reactions",
    "review-routes",
)

ROOT_FILES = (
    ".gitattributes",
    ".gitignore",
    "CITATION.cff",
    "CODE_OF_CONDUCT.md",
    "CONTRIBUTING.md",
    "LICENSE",
    "NOTICE",
    "README.md",
    "SECURITY.md",
    "THIRD_PARTY_NOTICES.md",
    "pyproject.toml",
    "requirements-dev.txt",
)

SKILL_FILES = (
    "SKILL.md",
    "agents/openai.yaml",
    "scripts/requirements.txt",
    "scripts/validate_output.py",
)

FORBIDDEN_NAMES = {
    ".DS_Store",
    ".env",
    "chemistry_skill_routing_f5.json",
    "chemistry_skill_routing_f5_agentplan_results.json",
    "reaction_audit_cases.json",
    "review_routes_expert_f5_candidates.json",
}

FORBIDDEN_PARTS = {
    ".pytest_cache",
    ".ruff_cache",
    ".venv",
    "__pycache__",
    "node_modules",
}

SECRET_PATTERNS = (
    re.compile(r"gh[pousr]_[A-Za-z0-9]{20,}"),
    re.compile(r"github_pat_[A-Za-z0-9_]{20,}"),
    re.compile(r"AKIA[0-9A-Z]{16}"),
    re.compile(r"sk-[A-Za-z0-9_-]{20,}"),
    re.compile(
        r"eyJ[A-Za-z0-9_-]{20,}\."
        r"[A-Za-z0-9_-]{20,}\."
        r"[A-Za-z0-9_-]{20,}"
    ),
)


def load_yaml(path: Path, errors: list[str]) -> object | None:
    try:
        return yaml.safe_load(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, yaml.YAMLError) as exc:
        errors.append(f"{path.relative_to(ROOT)}: invalid YAML: {exc}")
        return None


def frontmatter(path: Path, errors: list[str]) -> dict[str, object] | None:
    text = path.read_text(encoding="utf-8")
    if not text.startswith("---\n"):
        errors.append(f"{path.relative_to(ROOT)}: missing YAML frontmatter")
        return None
    try:
        raw = text.split("---\n", 2)[1]
    except IndexError:
        errors.append(f"{path.relative_to(ROOT)}: unclosed YAML frontmatter")
        return None
    try:
        value = yaml.safe_load(raw)
    except yaml.YAMLError as exc:
        errors.append(f"{path.relative_to(ROOT)}: invalid frontmatter: {exc}")
        return None
    if not isinstance(value, dict):
        errors.append(f"{path.relative_to(ROOT)}: frontmatter must be an object")
        return None
    return value


def pinned_requirements(path: Path, errors: list[str]) -> set[str]:
    packages: set[str] = set()
    for line_number, raw_line in enumerate(
        path.read_text(encoding="utf-8").splitlines(),
        start=1,
    ):
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        if not re.fullmatch(r"[A-Za-z0-9_.-]+==[A-Za-z0-9_.+-]+", line):
            errors.append(
                f"{path.relative_to(ROOT)}:{line_number}: "
                "dependency must use an exact == pin"
            )
            continue
        packages.add(line.lower())
    return packages


def validate_structure(errors: list[str]) -> None:
    for relative in ROOT_FILES:
        if not (ROOT / relative).is_file():
            errors.append(f"missing root file: {relative}")

    actual_skills = {
        path.name for path in SKILLS_ROOT.iterdir() if path.is_dir()
    }
    if actual_skills != set(SKILLS):
        errors.append(
            "skills directory mismatch: "
            f"expected={sorted(SKILLS)} actual={sorted(actual_skills)}"
        )

    for skill in SKILLS:
        skill_root = SKILLS_ROOT / skill
        for relative in SKILL_FILES:
            if not (skill_root / relative).is_file():
                errors.append(f"missing skill file: skills/{skill}/{relative}")
        scripts = list((skill_root / "scripts").glob("*.py"))
        if len(scripts) < 2:
            errors.append(f"skills/{skill}: expected processor and validator")


def validate_metadata(errors: list[str]) -> None:
    development_dependencies = pinned_requirements(
        ROOT / "requirements-dev.txt",
        errors,
    )
    for skill in SKILLS:
        skill_root = SKILLS_ROOT / skill
        metadata = frontmatter(skill_root / "SKILL.md", errors)
        if metadata is not None:
            if metadata.get("name") != skill:
                errors.append(f"skills/{skill}/SKILL.md: name must match directory")
            if set(metadata) != {"name", "description"}:
                errors.append(
                    f"skills/{skill}/SKILL.md: frontmatter permits only "
                    "name and description"
                )
        agent = load_yaml(skill_root / "agents" / "openai.yaml", errors)
        if not isinstance(agent, dict):
            errors.append(f"skills/{skill}/agents/openai.yaml: must be an object")
        skill_dependencies = pinned_requirements(
            skill_root / "scripts" / "requirements.txt",
            errors,
        )
        missing = skill_dependencies - development_dependencies
        if missing:
            errors.append(
                f"skills/{skill}: requirements-dev.txt is missing {sorted(missing)}"
            )

    citation = load_yaml(ROOT / "CITATION.cff", errors)
    if isinstance(citation, dict):
        required = {"cff-version", "message", "title", "authors"}
        if not required <= set(citation):
            errors.append("CITATION.cff: missing required CFF fields")
        if citation.get("license") != "Apache-2.0":
            errors.append("CITATION.cff: license must be Apache-2.0")

    for pattern in ("*.yaml", "*.yml"):
        for path in ROOT.rglob(pattern):
            if ".git" not in path.relative_to(ROOT).parts:
                load_yaml(path, errors)


def validate_local_links(errors: list[str]) -> None:
    link_pattern = re.compile(r"\[[^\]]+\]\(([^)]+)\)")
    for path in ROOT.rglob("*.md"):
        if ".git" in path.relative_to(ROOT).parts:
            continue
        text = path.read_text(encoding="utf-8")
        for raw_target in link_pattern.findall(text):
            target = raw_target.strip().split(maxsplit=1)[0].strip("<>")
            parsed = urlsplit(target)
            if parsed.scheme or target.startswith(("#", "mailto:")):
                continue
            relative_target = unquote(parsed.path)
            if not relative_target:
                continue
            resolved = (path.parent / relative_target).resolve()
            try:
                resolved.relative_to(ROOT)
            except ValueError:
                errors.append(
                    f"{path.relative_to(ROOT)}: link escapes repository: {target}"
                )
                continue
            if not resolved.exists():
                errors.append(
                    f"{path.relative_to(ROOT)}: broken local link: {target}"
                )


def validate_public_boundary(errors: list[str]) -> None:
    for path in ROOT.rglob("*"):
        relative = path.relative_to(ROOT)
        if any(part in FORBIDDEN_PARTS for part in relative.parts):
            continue
        if path.name in FORBIDDEN_NAMES:
            errors.append(f"forbidden private file: {relative}")
        if not path.is_file() or ".git" in relative.parts:
            continue
        if path.suffix.lower() in {".pyc", ".pyo"}:
            errors.append(f"forbidden compiled file: {relative}")
            continue
        if path.suffix.lower() not in {
            "",
            ".cff",
            ".md",
            ".py",
            ".toml",
            ".txt",
            ".yaml",
            ".yml",
        }:
            continue
        try:
            text = path.read_text(encoding="utf-8")
        except UnicodeError:
            errors.append(f"non-UTF-8 public text file: {relative}")
            continue
        user_path_marker = "/" + "Users" + "/"
        internal_path_marker = "byte" + "dance" + "/"
        if user_path_marker in text or internal_path_marker in text:
            errors.append(f"machine-specific or internal path in {relative}")
        for pattern in SECRET_PATTERNS:
            if pattern.search(text):
                errors.append(f"possible credential in {relative}")


def main() -> int:
    errors: list[str] = []
    validate_structure(errors)
    validate_metadata(errors)
    validate_local_links(errors)
    validate_public_boundary(errors)
    if errors:
        for error in errors:
            print(f"ERROR: {error}", file=sys.stderr)
        print(f"repository validation failed: {len(errors)} error(s)", file=sys.stderr)
        return 1
    print(f"repository validation passed: {len(SKILLS)} skills")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
