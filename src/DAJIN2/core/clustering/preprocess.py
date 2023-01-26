from __future__ import annotations
from collections import Counter
import numpy as np
from collections import defaultdict
from copy import deepcopy


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
        cnt = Counter(cssplit)
        if cnt["N"] == 0:  # No N
            continue
        if len(cnt) == 1 and cnt["N"] > 0:  # All N
            continue
        size = cnt["N"]
        del cnt["N"]
        samples = sampling(cnt, size)
        d_samples[i] = iter(samples)
    if not d_samples:
        return cssplits
    cssplits_replaced = deepcopy(cssplits)
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


def compress_insertion(cssplits: list[list[str]]) -> list[dict[str, int]]:
    """
    Insertion will be subdivided by mutations in the its sequence, so it is compressed as a '+I' to eliminate mutations.
    """
    cssplits_abstracted = []
    for cssplit in cssplits:
        for i, cs in enumerate(cssplit):
            if cs.startswith("+"):
                cssplit[i] = "+I" + cs.split("|")[-1]
        cssplits_abstracted.append(cssplit)
    return cssplits_abstracted

