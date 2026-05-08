"""Extract Python versions from pyproject.toml and write them to GITHUB_OUTPUT."""

from __future__ import annotations

import json
import os
import re
from pathlib import Path


def read_python_constraint(pyproject_path: str = "pyproject.toml") -> str:
    """Read the Python version constraint (e.g. '>=3.10,<3.13') from pyproject.toml."""
    text = Path(pyproject_path).read_text(encoding="utf-8")
    for pattern in [
        r'^\s*requires-python\s*=\s*"([^"]+)"',
        r'^\s*python\s*=\s*"([^"]+)"',
    ]:
        match = re.search(pattern, text, flags=re.MULTILINE)
        if match:
            return match.group(1).strip()
    raise RuntimeError("No Python version constraint found in pyproject.toml.")


def expand_python_versions(constraint: str) -> list[str]:
    """Expand a constraint like '>=3.10,<3.13' into ['3.10', '3.11', '3.12']."""
    operators = re.findall(r"(>=|<=|>|<|==|!=)\s*(\d+)\.(\d+)", constraint)
    if not operators:
        raise RuntimeError(f"Could not parse version constraint: {constraint!r}")

    # Collect all mentioned minor versions to determine the search range
    mentioned_minors = [int(minor) for _, _, minor in operators]
    min_minor = min(mentioned_minors)
    max_minor = max(mentioned_minors)

    major = int(operators[0][1])
    candidates = [f"{major}.{m}" for m in range(min_minor, max_minor + 1)]

    versions = []
    for version in candidates:
        v = tuple(map(int, version.split(".")))
        if all(_satisfies(v, op, int(maj), int(mn)) for op, maj, mn in operators):
            versions.append(version)
    return versions


def _satisfies(version: tuple[int, int], op: str, major: int, minor: int) -> bool:
    target = (major, minor)
    ops = {
        ">=": version >= target,
        ">": version > target,
        "<=": version <= target,
        "<": version < target,
        "==": version == target,
        "!=": version != target,
    }
    return ops[op]


def write_github_output(key: str, value: str) -> None:
    github_output = os.environ.get("GITHUB_OUTPUT")
    if github_output:
        with open(github_output, "a", encoding="utf-8") as fh:
            fh.write(f"{key}={value}\n")


def main() -> None:
    constraint = read_python_constraint()
    versions = expand_python_versions(constraint)
    versions_json = json.dumps(versions)
    print(f"Python versions: {versions_json}")
    write_github_output("versions", versions_json)
    write_github_output("latest_version", versions[-1])


if __name__ == "__main__":
    main()
