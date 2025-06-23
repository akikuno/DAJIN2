"""
Indel counting and normalization utilities shared between preprocess and consensus modules.
"""

from __future__ import annotations

from collections.abc import Iterator
from pathlib import Path

import numpy as np

from DAJIN2.utils import io


def count_indels(midsv_sample: Iterator[dict], sequence: str) -> dict[str, list[int]]:
    """Count indels from midsv data."""
    len_sequence = len(sequence)
    count = {"=": [0] * len_sequence, "+": [0] * len_sequence, "-": [0] * len_sequence, "*": [0] * len_sequence}
    for samp in midsv_sample:
        for i, cs in enumerate(samp["CSSPLIT"].split(",")):
            if cs == "N" or cs.islower():
                continue
            if cs.startswith("="):
                count["="][i] += 1
            elif cs.startswith("+"):
                count["+"][i] += 1
            elif cs.startswith("-"):
                count["-"][i] += 1
            elif cs.startswith("*"):
                count["*"][i] += 1
    return count


def normalize_indels(count: dict[str, list[int]]) -> dict[str, np.array]:
    """Normalize indel counts by total coverage."""
    count_normalized = {}
    match_count = np.array(count["="])
    for mut, indel_count in count.items():
        numerator = np.array(indel_count)
        denominator = numerator + match_count
        count_normalized[mut] = np.where(denominator == 0, 0, numerator / denominator * 100)
    return count_normalized


def minimize_mutation_counts(
    indels_control: dict[str, np.array], indels_sample: dict[str, np.array]
) -> dict[str, np.array]:
    """
    In cases where control has a larger value than sample, adjust the value of sample to match that of control.
    """
    indels_control_minimized = {}
    for mut in ["+", "-", "*"]:
        indels_control_minimized[mut] = np.minimum(indels_control[mut], indels_sample[mut])
    return indels_control_minimized


def summarize_indels(path_midsv: Path, sequence: str) -> tuple:
    """Returns indels, coverages, normalized indels, and kmer indels."""
    indels_count = count_indels(io.read_jsonl(path_midsv), sequence)
    indels_normalized = normalize_indels(indels_count)
    return indels_count, indels_normalized
