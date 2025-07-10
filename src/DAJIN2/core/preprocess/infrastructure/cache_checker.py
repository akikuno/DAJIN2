from __future__ import annotations

import hashlib
from pathlib import Path

from DAJIN2.utils import io


def exists_cached_hash(tempdir: Path, path: str) -> bool:
    path_cache_hash = Path(tempdir, "cache", "hash.txt")
    if path_cache_hash.exists():
        current_hash = hashlib.sha256(Path(path).read_bytes()).hexdigest()
        cashed_hash = path_cache_hash.read_text()
        if current_hash == cashed_hash:
            return True
    return False


def exists_cached_genome(tempdir: Path, genome: str) -> bool:
    path_cached_genome = Path(tempdir, "cache", "genome_coordinates.jsonl")
    if path_cached_genome.exists():
        cached_genome = list(io.read_jsonl(path_cached_genome))[0]
        if genome == cached_genome["genome"]:
            return True
    return False
