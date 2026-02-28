from __future__ import annotations

import hashlib
import json
from collections.abc import Iterable
from pathlib import Path


def to_allele_key(allele_name: str) -> str:
    """Return a stable MD5 hex key for an allele name used in internal paths."""
    return hashlib.md5(allele_name.encode("utf-8")).hexdigest()


def build_allele_name_map(allele_names: Iterable[str]) -> dict[str, dict[str, str]]:
    """
    Build bidirectional maps between original allele names and their MD5 keys to
    avoid 'OSError: [Errno 36] File name too long' due to long file paths.
    """
    name_to_key = {name: to_allele_key(name) for name in sorted(set(allele_names))}
    key_to_name = {key: name for name, key in name_to_key.items()}
    return {"name_to_key": name_to_key, "key_to_name": key_to_name}


def save_allele_name_map(tempdir: str | Path, allele_names: Iterable[str]) -> Path:
    """Persist the allele-name/key mapping JSON under tempdir/cache and return its path."""
    path_output = Path(tempdir, "cache", "allele_name_map.json")
    path_output.parent.mkdir(parents=True, exist_ok=True)
    path_output.write_text(json.dumps(build_allele_name_map(allele_names), indent=2), encoding="utf-8")
    return path_output
