from __future__ import annotations
from collections import Counter
import numpy as np
import re
from itertools import chain
from copy import deepcopy
from collections import defaultdict

from scipy import stats
from scipy.spatial.distance import cosine

from pathlib import Path

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


# def extract_indexes_with_both_ends_not_N(cssplits: str) -> tuple(int, int):
#     n_prefix = re.search(r"^(N,)+", cssplits)
#     left_idx = n_prefix.end() if n_prefix else 0
#     n_suffix = re.search(r"(,N)+$", cssplits)
#     right_idx = n_suffix.start() if n_suffix else len(cssplits)
#     # output index of splitted cssplits
#     left_idx = cssplits[:left_idx].count(",")
#     right_idx = cssplits.count(",") - cssplits[right_idx:].count(",")
#     return left_idx, right_idx


def call_count(cssplits: list[list[str]], indexes: list[tuple(int, int)]) -> dict[dict[str, int]]:
    """Count cssplits within 3-mer range.
    Args:
        cssplits (list[list[str]])
    Returns:
        list[dict[str, int]]: Both ends are counted as "N" to keep sequence length.
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


def call_percentage(cssplits: list[list[str]], counts: dict[dict[str, int]]) -> dict[dict[str:float]]:
    coverage = len(cssplits)
    percents = deepcopy(counts)
    for i, c in counts.items():
        for kmer, count in c.items():
            percents[i][kmer] = count / coverage * 100
    return percents


# call_percentage(counts)


# length = cssplits[0].count(",") + 1
# count_kmer = defaultdict(Counter)
# for cs in cssplits:
#     cs = cs.split(",")
#     for i in range(1, length - 1):
#         print(cs[i - 1], cs[i], cs[i + 1])
#         kmer = ",".join([cs[i - 1], cs[i], cs[i + 1]])
#         count_kmer[i] += Counter([kmer])


# # cssplits = "N,N,N,N,=A,=C,=G,=T"
# cssplits = "N,N,=A,=C,=G,=T,N,N"
# left_idx, right_idx = extract_indexes_with_both_ends_not_N(cssplits)
# # cssplits.split(",")[left_idx]
# # cssplits.split(",")[right_idx]

# cssplits_sample = ["=A,=C,=G,=T", "*AT,-C,+A|G,=t"]

# a = [0.01, 0.02, 0.03]
# b = [0.97, 0.98, 0.99]
# 1 - cosine(a, b)

# _, pval = stats.ttest_ind(a, b, equal_var=False)
# pval
