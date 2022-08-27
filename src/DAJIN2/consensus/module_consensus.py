from __future__ import annotations
from collections import Counter


def call(cssplits: list[str]) -> list[str]:
    cs = [cs.split(",") for cs in cssplits]
    cons = []
    for c in zip(*cs):
        cons.append(Counter(c).most_common()[0][0])
    return ",".join(cons)
