from __future__ import annotations

import json
from pathlib import Path

from rapidfuzz import process
from rapidfuzz.distance import DamerauLevenshtein

from DAJIN2.utils import io


def extract_unique_sv(
    fasta_sv_alleles: dict[str, str], FASTA_ALLELES: dict[str, str], base_num: int = 10
) -> dict[str, str]:
    """
    Extract unique SVs (insertions/inversions) alleles if they are dissimilar to the FASTA_ALLELES input by the user.
    "Unique SV alleles" are defined as sequences that have a difference of more than 10 bases compared to the sequences in FASTA_ALLELES
    """
    fasta_sv_alleles_unique = fasta_sv_alleles.copy()

    # Remove SV alleles that are similar to the FASTA_ALLELES input by the user
    to_delete = set()
    for key, seq in fasta_sv_alleles_unique.items():
        _, distances, _ = zip(*process.extract_iter(seq, FASTA_ALLELES.values(), scorer=DamerauLevenshtein.distance))
        if any(d < base_num for d in distances):
            to_delete.add(key)

    # Remove SV alleles that are similar to each other
    for key, seq in fasta_sv_alleles_unique.items():
        if key in to_delete:
            continue
        _, distances, _ = zip(
            *process.extract_iter(seq, fasta_sv_alleles_unique.values(), scorer=DamerauLevenshtein.distance)
        )
        similar_index = {
            k if d < base_num else None for k, d in zip(fasta_sv_alleles_unique.keys(), distances) if k != key
        }
        to_delete |= similar_index

    return {k: v for k, v in fasta_sv_alleles_unique.items() if k not in to_delete}


###########################################################
# Update keys to avoid duplicating user-specified alleles
###########################################################


def _check_duplicates_of_sets(set1: set[str], set2: set[str]) -> bool:
    """if True, there are duplicates."""
    union_set = set1.union(set2)
    total_elements = len(set1) + len(set2)
    return len(union_set) != total_elements


def add_unique_allele_keys(
    fasta_sv_alleles: dict[str, str],
    FASTA_ALLELES: dict[str, set],
    key: str,
    internal_suffix: str = "",
    display_prefix: str = "DAJIN_",
) -> tuple[dict[str, str], dict[str, str]]:
    """
    Update keys to avoid duplicating user-specified alleles.
    If the allele 'insertion01' exists in FASTA_ALLELES, increment the digits.
    (insertion01 -> insertion001 -> insertion0001...)

    Returns:
        tuple[dict[str, str], dict[str, str]]:
            - A dictionary whose keys are internal allele names (optionally suffixed with `internal_suffix`)
            - A mapping from the internal allele name to the display allele name (prefixed with `display_prefix`)
    """
    user_defined_alleles = set(FASTA_ALLELES)
    key_duplicated_alleles = {allele for allele in user_defined_alleles if key in allele}

    if key_duplicated_alleles == set():
        base_names = [f"{key}{(i + 1):02}" for i, _ in enumerate(fasta_sv_alleles.values())]
    else:
        num_digits = 3  # 001
        while True:
            base_names = [f"{key}{(i + 1):0{num_digits}}" for i, _ in enumerate(fasta_sv_alleles)]
            if not _check_duplicates_of_sets(set(base_names), key_duplicated_alleles):
                break
            num_digits += 1

    internal_names = [f"{name}__{internal_suffix}" if internal_suffix else name for name in base_names]
    display_names = [f"{display_prefix}{name}" if display_prefix else name for name in base_names]

    return dict(zip(internal_names, fasta_sv_alleles.values())), dict(zip(internal_names, display_names))


###########################################################
# Save fasta and midsv files
###########################################################


def save_fasta(TEMPDIR: Path | str, SAMPLE_NAME: str, fasta_sv_alleles: dict[str, str]) -> None:
    Path(TEMPDIR, SAMPLE_NAME, "fasta").mkdir(parents=True, exist_ok=True)
    for header, seq in fasta_sv_alleles.items():
        Path(TEMPDIR, SAMPLE_NAME, "fasta", f"{header}.fasta").write_text(f">{header}\n{seq}\n")


def save_midsv(TEMPDIR: Path | str, SAMPLE_NAME: str, midsv_sv_alleles: dict[str, list[str]]) -> None:
    Path(TEMPDIR, SAMPLE_NAME, "midsv").mkdir(parents=True, exist_ok=True)
    for header, midsv_tag in midsv_sv_alleles.items():
        io.write_jsonl(midsv_tag, Path(TEMPDIR, SAMPLE_NAME, "midsv", f"consensus_{header}.jsonl"))


def save_sv_name_map(TEMPDIR: Path | str, SAMPLE_NAME: str, name_map: dict[str, str]) -> None:
    """Persist internal->display SV allele name mapping for later report generation."""
    path_map = Path(TEMPDIR, SAMPLE_NAME, "fasta", "sv_name_map.json")
    path_map.parent.mkdir(parents=True, exist_ok=True)
    path_map.write_text(json.dumps(name_map, indent=2))


def load_sv_name_map(TEMPDIR: Path | str, SAMPLE_NAME: str) -> dict[str, str]:
    """Load internal->display SV allele name mapping if it exists."""
    path_map = Path(TEMPDIR, SAMPLE_NAME, "fasta", "sv_name_map.json")
    if not path_map.exists():
        return {}
    return json.loads(path_map.read_text())


def invert_sv_name_map(name_map: dict[str, str]) -> dict[str, str]:
    """Invert {internal: display} to {display: internal}."""
    return {display: internal for internal, display in name_map.items()}
