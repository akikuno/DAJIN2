from __future__ import annotations

import re

import numpy as np


def get_repeat_regions(sequence: str, loci: set[int]) -> list[tuple[int, int]]:
    """
    Find homopolymers in the sequence but discard them that
    are adjacent to candidate mutation loci because they are
    likely to be covered by the real mutations
    """
    pattern = r"A{4,}|C{4,}|G{4,}|T{4,}|N{4,}"
    repeat_regions = []
    for start, end in (match.span() for match in re.finditer(pattern, sequence)):
        if not (start - 1 in loci and end + 1 in loci):
            repeat_regions.append((start, end))
    return repeat_regions


def cosine_simirarity(X, Y) -> float:
    return np.dot(X, Y) / (np.linalg.norm(X) * np.linalg.norm(Y))


###########################################################
# main
###########################################################


def extract_errors(
    sequence: str,
    indels_normalized_sample: dict[str, np.array],
    indels_normalized_control: dict[str, np.array],
    candidate_loci: dict[set],
) -> dict[str, set(int)]:
    errors_in_homopolymer = dict()
    for mut in ["+", "-", "*"]:
        repeat_regions = get_repeat_regions(sequence, candidate_loci[mut])
        if len(repeat_regions) == 0:
            errors_in_homopolymer[mut] = set()
            continue
        errors = set()
        for start, end in repeat_regions:
            x = indels_normalized_sample[mut][start:end]
            y = indels_normalized_control[mut][start:end]
            if cosine_simirarity(x, y) > 0.95:
                errors.update(range(start, end + 1))
        errors_in_homopolymer[mut] = errors

    return errors_in_homopolymer
