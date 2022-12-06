from __future__ import annotations
from collections import Counter
import numpy as np
from collections import defaultdict


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


def replace_both_ends_n_by_sampling(transpose_cssplits: list[list[str]]):
    d_samples = defaultdict(iter)
    for i, cssplits in enumerate(transpose_cssplits):
        cnt = Counter(cssplits)
        if cnt["N"] == 0:  # No N
            continue
        if len(cnt) == 1 and cnt["N"] > 0:  # All N
            continue
        size = cnt["N"]
        del cnt["N"]
        samples = sampling(cnt, size)
        d_samples[i] = iter(samples)
    cssplits_replaced = [list(cs) for cs in zip(*transpose_cssplits)]
    for i, cssplits in enumerate(cssplits_replaced):
        for j, cs in enumerate(cssplits):
            if cs != "N":
                break
            cssplits_replaced[i][j] = next(d_samples[j])
    for i, cssplits in enumerate(cssplits_replaced):
        cssplits = cssplits[::-1]
        for j, cs in enumerate(cssplits):
            if cs != "N":
                break
            try:
                cssplits_replaced[i][len(cssplits) - 1 - j] = next(d_samples[len(cssplits) - 1 - j])
            except TypeError:
                pass
    cssplits_replaced = cssplits_replaced[::-1]
    cssplits_replaced = [list(cs) for cs in zip(*cssplits_replaced)]
    return cssplits_replaced


###############################################################################
# main
###############################################################################


def replace_both_ends_n(cssplits_control: list[list[str]]):
    transpose_control = transpose(cssplits_control)
    # Make count matrix
    count_control = call_count(transpose_control)
    # Distribute N
    transpose_compensated_control = replace_both_ends_n_by_sampling(count_control)
    return transpose(transpose_compensated_control)


import midsv

sam = midsv.read_sam("tmp_ins.sam")
x = midsv.transform(sam, midsv=False, qscore=False)

x[0]["CSSPLIT"].count(",")
