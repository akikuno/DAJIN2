from __future__ import annotations

import re

import numpy as np
from numpy.linalg import norm


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


def cosine_similarity_safe(x: np.ndarray, y: np.ndarray) -> float:
    """Compute cosine similarity while avoiding divide-by-zero warnings."""
    if x.size == 0 or y.size == 0:
        return 0.0

    denominator = norm(x) * norm(y)
    if denominator == 0 or np.isnan(denominator):
        return 0.0

    return float(np.dot(x, y) / denominator)


###########################################################
# main
###########################################################


def extract_sequence_errors_in_homopolymer_loci(
    sequence: str,
    indels_normalized_sample: dict[str, np.array],
    indels_normalized_control: dict[str, np.array],
    anomal_loci: dict[set],
) -> dict[str, set[int]]:
    sequence_errors_in_homopolymer = {}
    for mut in ["+", "-", "*"]:
        repeat_regions = get_repeat_regions(sequence, anomal_loci[mut])
        if len(repeat_regions) == 0:
            sequence_errors_in_homopolymer[mut] = set()
            continue
        sequence_errors = set()
        for start, end in repeat_regions:
            x = np.array(indels_normalized_sample[mut][start:end], dtype=float)
            y = np.array(indels_normalized_control[mut][start:end], dtype=float)

            # Treat variants >= 3% in the control as homopolymer sequence errors and replace with sample values.
            # Replacing with sample values increases the cosine similarity.
            mask = y >= 3
            y[mask] = x[mask]

            cos_sim = cosine_similarity_safe(x, y)
            # spearman_corr, _ = spearmanr(x, y)
            if cos_sim > 0.8:
                sequence_errors.update(range(start, end + 1))

        sequence_errors_in_homopolymer[mut] = sequence_errors

    return sequence_errors_in_homopolymer
