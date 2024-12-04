from __future__ import annotations

from pathlib import Path

from rapidfuzz import process
from rapidfuzz.distance import DamerauLevenshtein


def extract_unique_sv(fasta_sv_alleles: dict[str, str], FASTA_ALLELES: dict[str, str]) -> dict[str, str]:
    """
    Extract unique SVs (insertions/inversions) alleles if they are dissimilar to the FASTA_ALLELES input by the user.
    "Unique SV alleles" are defined as sequences that have a difference of more than 10 bases compared to the sequences in FASTA_ALLELES
    """
    fasta_sv_alleles_unique = fasta_sv_alleles.copy()

    # Remove SV alleles that are similar to the FASTA_ALLELES input by the user
    to_delete = set()
    for key, seq in fasta_sv_alleles_unique.items():
        _, distances, _ = zip(*process.extract_iter(seq, FASTA_ALLELES.values(), scorer=DamerauLevenshtein.distance))
        if any(d < 10 for d in distances):
            to_delete.add(key)

    # Remove SV alleles that are similar to each other
    for key, seq in fasta_sv_alleles_unique.items():
        if key in to_delete:
            continue
        _, distances, _ = zip(
            *process.extract_iter(seq, fasta_sv_alleles_unique.values(), scorer=DamerauLevenshtein.distance)
        )
        similar_index = {k if d < 10 else None for k, d in zip(fasta_sv_alleles_unique.keys(), distances) if k != key}
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


def add_unique_allele_keys(fasta_sv_alleles: dict[str, str], FASTA_ALLELES: dict[str, set], key: str) -> dict[str, str]:
    """
    Update keys to avoid duplicating user-specified alleles.
    If the allele 'insertion1' exists in FASTA_ALLELES, increment the digits.
    (insertion1 -> insertion01 -> insertion001...)
    """
    user_defined_alleles = set(FASTA_ALLELES)
    key_duplicated_alleles = {allele for allele in user_defined_alleles if key in allele}
    if key_duplicated_alleles == set():
        return {f"{key}{i+1}": value for i, value in enumerate(fasta_sv_alleles.values())}

    key_candidate_alleles = set()
    digits = 2
    while _check_duplicates_of_sets(key_candidate_alleles, key_duplicated_alleles):
        key_candidate_alleles = {f"{key}{i+1:0{digits}}" for i, _ in enumerate(fasta_sv_alleles)}
        digits += 1

    return dict(zip(key_candidate_alleles, fasta_sv_alleles.values()))


###########################################################
# Save cstag and fasta
###########################################################


def save_fasta(TEMPDIR: Path | str, SAMPLE_NAME: str, fasta_sv_alleles: dict[str, str]) -> None:
    for header, seq in fasta_sv_alleles.items():
        Path(TEMPDIR, SAMPLE_NAME, "fasta", f"{header}.fasta").write_text(f">{header}\n{seq}\n")


def save_cstag(TEMPDIR: Path | str, SAMPLE_NAME: str, cstag_sv_alleles: dict[str, str]) -> None:
    for header, cs_tag in cstag_sv_alleles.items():
        Path(TEMPDIR, SAMPLE_NAME, "cstag", f"{header}.txt").write_text(cs_tag + "\n")
