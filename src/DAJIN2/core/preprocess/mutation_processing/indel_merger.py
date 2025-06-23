"""
Utilities for merging contiguous indels and handling mutation loci.
"""

from __future__ import annotations

import bisect
from collections import defaultdict

import numpy as np


def count_elements_within_range(arr, lower_bound, upper_bound):
    """
    Counts the number of elements within a specified range in a sorted array.
    """
    start_index = bisect.bisect_left(arr, lower_bound)
    end_index = bisect.bisect_right(arr, upper_bound)
    return end_index - start_index


def merge_index_of_consecutive_indel(mutation_loci: dict[str, set[int]]) -> dict[str, set[int]]:
    """Treat as contiguous indels if there are insertions/deletions within five bases of each other"""
    mutation_loci_merged = {}

    # Reflect point mutations as they are
    mutation_loci_merged["*"] = mutation_loci["*"]

    # Merge if indels are within 10 bases
    for mut in ["+", "-"]:
        idx_indel = sorted(mutation_loci[mut])
        idx_indel_merged = set(idx_indel)
        for i in range(len(idx_indel) - 1):
            idx_1 = idx_indel[i]
            idx_2 = idx_indel[i + 1]
            # If everything from idx_1 to idx_2 is already considered as indels, then skip it.
            if count_elements_within_range(idx_indel, idx_1 + 1, idx_2 - 1) == idx_2 - idx_1 + 1:
                continue
            if idx_1 + 10 > idx_2:
                for i in range(idx_1 + 1, idx_2):
                    idx_indel_merged.add(i)
        mutation_loci_merged[mut] = idx_indel_merged

    # Additional logic for mutation enrichment within 10 bases on both ends
    for mut in ["+", "-"]:
        idx_indel = sorted(mutation_loci_merged[mut])
        idx_indel_merged = set(idx_indel)
        for i in range(len(idx_indel) - 1):
            idx_1 = idx_indel[i]
            idx_2 = idx_indel[i + 1]
            # If everything from idx_1 to idx_2 is already considered as indels, then skip it.
            if count_elements_within_range(idx_indel, idx_1 + 1, idx_2 - 1) == idx_2 - idx_1 + 1:
                continue
            # If the distance between idx_1 and idx_2 is more than 20 bases, then skip it.
            if idx_1 + 20 < idx_2:
                continue
            count_left = count_elements_within_range(idx_indel, idx_1 - 11, idx_1 - 1)
            count_right = count_elements_within_range(idx_indel, idx_2 + 1, idx_2 + 11)

            # If 8 out of the 10 bases at both ends are indels,
            # then everything from idx_1 to idx_2 will be considered as indels.

            if count_left >= 8 and count_right >= 8:
                for i in range(idx_1 + 1, idx_2):
                    idx_indel_merged.add(i)
        mutation_loci_merged[mut] = idx_indel_merged

    return mutation_loci_merged


def split_kmer(indels: dict[str, np.array], kmer: int = 11) -> dict[str, np.array]:
    """Split indel data into k-mer windows."""
    results = defaultdict(list)
    center = kmer // 2
    for mut, value in indels.items():
        for i in range(len(value)):
            if center <= i <= len(value) - center:
                start = i - center
                if kmer % 2 == 0:
                    end = i + center
                else:
                    end = i + center + 1
                results[mut].append(value[start:end])
            else:
                results[mut].append(np.array([0] * kmer))
    return results


def discard_errors(loci: dict[str, set[int]], errors: dict[str, set[int]]) -> dict[str, set[int]]:
    """Remove detected errors from the candidate mutation loci."""
    return {mut: loci[mut] - errors[mut] for mut in {"+", "-", "*"}}


def merge_loci(dissimilar_loci: dict[str, set], anomal_loci: dict[str, set]) -> dict[str, set]:
    """Merge dissimilar loci and anomalous loci."""
    mutation_loci = {}
    for mut in {"+", "-", "*"}:
        mutation_loci[mut] = dissimilar_loci[mut] | anomal_loci[mut]
    return mutation_loci


def add_knockin_loci(candidate_loci: dict[str, set], knockin_loci: set):
    """Add knockin loci to candidate mutation loci."""
    mutation_loci = {}
    for mut in {"+", "-", "*"}:
        mutation_loci[mut] = candidate_loci[mut] | knockin_loci
    return mutation_loci


def transpose_mutation_loci(mutation_loci: dict[str, set[int]], sequence: str) -> list[set[str]]:
    """Transpose mutation loci from mutation-type indexed to position-indexed format."""
    len_sequence = len(sequence)
    mutation_loci_transposed = [set() for _ in range(len_sequence)]
    for mut, idx_mutation in mutation_loci.items():
        for i, loci in enumerate(mutation_loci_transposed):
            if i in idx_mutation:
                loci.add(mut)
    return mutation_loci_transposed
