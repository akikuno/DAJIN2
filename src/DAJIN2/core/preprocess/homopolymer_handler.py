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


def cosine_similarity(X, Y) -> float:
    # Add 1e-6 to avoid division by zero when calculating cosine similarity
    X += 1e-6
    Y += 1e-6
    return float(np.dot(X, Y) / (np.linalg.norm(X) * np.linalg.norm(Y)))


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
            x = np.array(indels_normalized_sample[mut][start:end])
            y = np.array(indels_normalized_control[mut][start:end])

            # Scaling data to [0, 1] for cosine similarity

            # Check if the range of x is zero
            if x.max() - x.min() == 0:
                x_scaled = np.zeros_like(x)
            else:
                x_scaled = (x - x.min()) / (x.max() - x.min())

            # Check if the range of y is zero
            if y.max() - y.min() == 0:
                y_scaled = np.zeros_like(y)
            else:
                y_scaled = (y - y.min()) / (y.max() - y.min())

            if cosine_similarity(x_scaled, y_scaled) > 0.95:
                sequence_errors.update(range(start, end + 1))

        sequence_errors_in_homopolymer[mut] = sequence_errors

    return sequence_errors_in_homopolymer
