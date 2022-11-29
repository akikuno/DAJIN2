from __future__ import annotations
from copy import deepcopy
from collections import defaultdict


def reorder_labels(labels: list[int], start: int = 0) -> list[int]:
    labels_ordered = deepcopy(labels)
    num = start
    d = defaultdict(int)
    for i, l in enumerate(labels_ordered):
        if not d[l]:
            num += 1
            d[l] = num
        labels_ordered[i] = d[l]
    return labels_ordered

