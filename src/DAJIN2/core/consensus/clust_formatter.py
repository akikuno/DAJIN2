from __future__ import annotations

import random
from itertools import groupby
from collections import defaultdict


def remove_minor_alleles(clust_sample: list[dict]) -> list[dict]:
    """Remove minor alleles from the clustering sample."""
    min_sample_size = max(5, int(len(clust_sample) * 0.5 // 100))

    counts = defaultdict(int)
    for c in clust_sample:
        counts[c["LABEL"]] += 1

    minor_labels = {label for label, count in counts.items() if count < min_sample_size}

    return [c for c in clust_sample if c["LABEL"] not in minor_labels]


def downsample_by_label(clust_sample: list[dict], num: int = 1000) -> list[dict]:
    random.seed(1)

    clust_downsampled = []
    clust_sample.sort(key=lambda x: x["LABEL"])
    for _, group in groupby(clust_sample, key=lambda x: x["LABEL"]):
        samples = list(group)
        num_samples = min(len(samples), num)
        clust_downsampled.extend(random.sample(samples, num_samples))

    return clust_downsampled
