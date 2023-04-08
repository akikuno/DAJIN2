from __future__ import annotations

from collections import Counter, defaultdict

import numpy as np


def transpose(cssplits):
    return [list(cs) for cs in zip(*cssplits)]


def call_count(transpose_cssplits: list[list[str]]) -> list[dict[str:int]]:
    cssplit_counts = []
    for cssplit in transpose_cssplits:
        count = Counter(cssplit)
        count = dict(count)
        cssplit_counts.append(count)
    return cssplit_counts


def sampling(cnt: Counter, size: int):
    elements = []
    probs = []
    coverage = sum(cnt.values())
    for key, val in cnt.items():
        elements.append(key)
        probs.append(val / coverage)
    np.random.seed(1)
    samples = np.random.choice(a=elements, size=size, p=probs)
    return samples


###############################################################################
# main
###############################################################################


def replace_both_ends_n(cssplits: list[list[str]]):
    transpose_cssplits = [list(cs) for cs in zip(*cssplits)]
    d_samples = defaultdict(iter)
    for i, cssplit in enumerate(transpose_cssplits):
        flag_no_N = all(True if cs != "N" else False for cs in cssplit)
        if flag_no_N:
                continue
        flag_all_N = all(True if cs == "N" else False for cs in cssplit)
        cnt = Counter(cssplit)
        size = cnt["N"]
        if not flag_all_N:
            del cnt["N"]
        samples = sampling(cnt, size)
        d_samples[i] = iter(samples)
    cssplits_replaced = cssplits.copy()
    for i, cssplit in enumerate(cssplits_replaced):
        for j, cs in enumerate(cssplit):
            if cs != "N":
                break
            cssplits_replaced[i][j] = next(d_samples[j])
    for i, cssplit in enumerate(cssplits_replaced):
        cssplit = cssplit[::-1]
        for j, cs in enumerate(cssplit):
            if cs != "N":
                break
            try:
                cssplits_replaced[i][len(cssplit) - 1 - j] = next(d_samples[len(cssplit) - 1 - j])
            except TypeError:
                pass
    return cssplits_replaced


