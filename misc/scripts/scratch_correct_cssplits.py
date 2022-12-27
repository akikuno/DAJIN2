from __future__ import annotations
from collections import Counter
import re
from itertools import chain
from copy import deepcopy
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy import stats
from scipy.spatial.distance import cosine

import midsv

"""
- 5-merでコントロールと比較して、コントロールにもあるプロファイルの場合はマッチで補正する
    - 断片のリードで、端から続くNは無視
    - 補正する前にすでにすべてがマッチなら即continue
"""


def extract_indexes_with_both_ends_not_N(cssplits: list[list[str]]) -> list[tuple[int, int]]:
    indexes = []
    for cssplit in cssplits:
        cssplit = ",".join(cssplit)
        n_prefix = re.search(r"^(N,)+", cssplit)
        left_idx = n_prefix.end() if n_prefix else 0
        n_suffix = re.search(r"(,N)+$", cssplit)
        right_idx = n_suffix.start() if n_suffix else len(cssplit)
        # output index of splitted cssplit
        left_idx = cssplit[:left_idx].count(",")
        right_idx = cssplit.count(",") - cssplit[right_idx:].count(",")
        indexes.append((left_idx, right_idx))
    return indexes


def call_count(cssplits: list[list[str]], indexes: list[tuple(int, int)]) -> dict[dict[str, int]]:
    """Count cssplits within 3-mer range.
    """
    count_kmer = defaultdict(Counter)
    for cssplit, idx in zip(cssplits, indexes):
        left_idx, right_idx = idx
        for i in range(left_idx + 1, right_idx):
            if cssplit[i].startswith("="):
                continue
            kmer = ",".join([cssplit[i - 1], cssplit[i], cssplit[i + 1]])
            count_kmer[i] += Counter([kmer])
    counts = {i: dict(count_kmer[i]) for i in count_kmer.keys()}
    return counts


def call_percentage(cssplits: list[list[str]], counts: dict[dict[str, int]]) -> dict[dict[str, float]]:
    coverage = len(cssplits)
    percents = deepcopy(counts)
    for i, c in counts.items():
        for kmer, count in c.items():
            percents[i][kmer] = count / coverage * 100
    return percents


def subtract_percentage(percent_sample: dict, percent_control: dict) -> dict[dict[str, float]]:
    percent_subtracted = deepcopy(percent_sample)
    for i, samp in percent_sample.items():
        cont = percent_control.get(i)
        for kmer_samp, per_samp in samp.items():
            if cont.get(kmer_samp):
                percent_subtracted[i][kmer_samp] = per_samp - cont[kmer_samp]
    return percent_subtracted


def select_candidate_mutation(percent_subtracted: dict, threshold: float = 0.5) -> dict[dict[str, float]]:
    candidate_mutation = dict()
    for i, samp in percent_subtracted.items():
        mutation = {k for k, v in samp.items() if v > threshold}
        candidate_mutation.update({i: mutation})
    return candidate_mutation


def update_cssplits(cssplits: list, sequence: str, candidate_mutation: dict) -> list[list[str]]:
    cssplits_update = deepcopy(cssplits)
    for j, mutation in candidate_mutation.items():
        for i, cssplit in enumerate(cssplits):
            kmer = ",".join([cssplit[j - 1], cssplit[j], cssplit[j + 1]])
            if kmer in mutation:
                continue
            else:
                cssplits_update[i][j] = "=" + sequence[j]
    return cssplits_update

