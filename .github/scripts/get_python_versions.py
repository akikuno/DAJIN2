from __future__ import annotations

import json
import os
import re
import subprocess
import sys
import urllib.request
from pathlib import Path

PYTHON_DEVGUIDE_URL = "https://devguide.python.org/versions/"
PACKAGES = ("mappy", "pysam")
SUBDIRS = ("linux-64", "osx-64")
DEFAULT_VERSIONS = ["3.10", "3.11", "3.12", "3.13"]


def sort_versions(versions: list[str] | set[str]) -> list[str]:
    return sorted(set(versions), key=lambda v: tuple(map(int, v.split("."))))


def normalize_py_tag(tag: str) -> str:
    """
    Convert:
      py39   -> 3.9
      py310  -> 3.10
      py311  -> 3.11
    """
    digits = tag.replace("py", "")
    if len(digits) == 2:
        return f"{digits[0]}.{digits[1]}"
    if len(digits) == 3:
        return f"{digits[0]}.{digits[1:]}"
    raise ValueError(f"Unexpected Python build tag: {tag}")


def fetch_supported_python_versions() -> list[str]:
    with urllib.request.urlopen(PYTHON_DEVGUIDE_URL, timeout=30) as response:
        html = response.read().decode("utf-8", errors="replace")

    match = re.search(
        r"Supported versions(.*?)Unsupported versions",
        html,
        flags=re.DOTALL,
    )
    if not match:
        raise RuntimeError("Could not find the 'Supported versions' section in the Python devguide.")

    section = match.group(1)
    versions = re.findall(r">(\d+\.\d+)<", section)
    versions = sort_versions(versions)

    if not versions:
        raise RuntimeError("No supported Python versions were extracted from the Python devguide.")

    return versions


def read_python_constraint() -> str:
    pyproject_text = Path("pyproject.toml").read_text(encoding="utf-8")

    try:
        import tomllib
    except ModuleNotFoundError:
        tomllib = None

    if tomllib is not None:
        data = tomllib.loads(pyproject_text)

        requires_python = data.get("project", {}).get("requires-python")
        if isinstance(requires_python, str) and requires_python.strip():
            return requires_python.strip()

        poetry_python = data.get("tool", {}).get("poetry", {}).get("dependencies", {}).get("python")
        if isinstance(poetry_python, str) and poetry_python.strip():
            return poetry_python.strip()

    patterns = [
        r'^\s*requires-python\s*=\s*"([^"]+)"',
        r'^\s*python\s*=\s*"([^"]+)"',
    ]
    for pattern in patterns:
        match = re.search(pattern, pyproject_text, flags=re.MULTILINE)
        if match:
            return match.group(1).strip()

    raise RuntimeError("No Python version constraint was found in pyproject.toml.")


def version_satisfies_simple(version: str, constraint: str) -> bool:
    current = tuple(map(int, version.split(".")))
    parts = re.findall(r"(<=|>=|==|!=|<|>)\s*(\d+(?:\.\d+)*)", constraint)

    if not parts:
        return True

    for operator, raw_target in parts:
        target = tuple(map(int, raw_target.split(".")))

        if operator == ">=" and not (current >= target):
            return False
        if operator == ">" and not (current > target):
            return False
        if operator == "<=" and not (current <= target):
            return False
        if operator == "<" and not (current < target):
            return False
        if operator == "==" and not (current == target):
            return False
        if operator == "!=" and not (current != target):
            return False

    return True


def filter_versions_by_constraint(versions: list[str], constraint: str) -> list[str]:
    try:
        from packaging.specifiers import SpecifierSet
        from packaging.version import Version

        specifier = SpecifierSet(constraint)
        return [v for v in versions if Version(v) in specifier]
    except Exception:
        return [v for v in versions if version_satisfies_simple(v, constraint)]


def conda_search_records(package: str, subdir: str) -> list[dict]:
    command = [
        "conda",
        "search",
        "--json",
        "--override-channels",
        "-c",
        "bioconda",
        "--subdir",
        subdir,
        package,
    ]
    result = subprocess.run(
        command,
        check=True,
        capture_output=True,
        text=True,
    )
    payload = json.loads(result.stdout)

    records = payload.get(package, [])
    if not records:
        raise RuntimeError(f"No conda records were found for package={package!r}, subdir={subdir!r}.")

    return records


def extract_python_versions_from_records(records: list[dict]) -> list[str]:
    versions = set()

    for record in records:
        build = record.get("build", "")
        match = re.search(r"py(\d{2,3})", build)
        if match:
            versions.add(normalize_py_tag(match.group(0)))

    versions = sort_versions(versions)

    if not versions:
        raise RuntimeError("No Python build tags were extracted from conda records.")

    return versions


def get_package_supported_versions(package: str, subdir: str) -> set[str]:
    records = conda_search_records(package, subdir)
    return set(extract_python_versions_from_records(records))


def get_shared_conda_versions(packages: tuple[str, ...], subdirs: tuple[str, ...]) -> list[str]:
    shared: set[str] | None = None

    for subdir in subdirs:
        for package in packages:
            supported = get_package_supported_versions(package, subdir)
            shared = supported if shared is None else (shared & supported)

    return sort_versions(shared or [])


def write_github_output(key: str, value: str) -> None:
    github_output = os.environ.get("GITHUB_OUTPUT")
    if github_output:
        with open(github_output, "a", encoding="utf-8") as fh:
            fh.write(f"{key}={value}\n")


def main() -> int:
    try:
        active_versions = fetch_supported_python_versions()
        python_constraint = read_python_constraint()
        project_versions = filter_versions_by_constraint(active_versions, python_constraint)
        conda_versions = get_shared_conda_versions(PACKAGES, SUBDIRS)

        final_versions = sort_versions(set(project_versions) & set(conda_versions))

        if not final_versions:
            raise RuntimeError(
                "No Python versions matched all filters. "
                f"active_versions={active_versions}, "
                f"python_constraint={python_constraint!r}, "
                f"project_versions={project_versions}, "
                f"conda_versions={conda_versions}"
            )

    except Exception as exc:
        print(f"::warning::{exc}", file=sys.stderr)
        final_versions = DEFAULT_VERSIONS

    versions_json = json.dumps(final_versions)
    print(f"Selected Python versions: {versions_json}")
    write_github_output("versions", versions_json)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
