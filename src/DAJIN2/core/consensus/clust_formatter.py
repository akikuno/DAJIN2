from __future__ import annotations

from itertools import groupby


def subset_clust(clust_sample: list[dict], num: int = 1000) -> list[dict]:
    clust_subset_sample = []
    clust_sample.sort(key=lambda x: x["LABEL"])
    for _, group in groupby(clust_sample, key=lambda x: x["LABEL"]):
        clust_subset_sample.extend(list(group)[:num])
    return clust_subset_sample
