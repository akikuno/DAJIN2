from __future__ import annotations
from copy import deepcopy
from collections import defaultdict, Counter


def merge_mixed_cluster(labels_control: list[int], labels_sample: list[int]) -> list[int]:
    labels_all = labels_sample + labels_control
    labels_merged = deepcopy(labels_sample)
    coverage_control = len(labels_control)
    labels_percent_control = [0] * len(Counter(labels_all).keys())
    for i, label in enumerate(labels_control):
        labels_percent_control[label] += 1 / coverage_control * 100
    labels_mixed = {i for i, per in enumerate(labels_percent_control) if per > 0.5}
    max_label = max(labels_merged) + 1
    for i, label in enumerate(labels_sample):
        if label in labels_mixed:
            labels_merged[i] = max_label
    return labels_merged


def merge_minor_cluster(labels_sample: list[int]) -> list[int]:
    total_num = sum(Counter(labels_sample).values())
    cnt = Counter(labels_sample)
    labels_minor = {label for label, num in cnt.items() if num / total_num * 100 < 0.5}
    max_label = max(labels_sample) + 1
    labels_merged = deepcopy(labels_sample)
    for i, label in enumerate(labels_sample):
        if label in labels_minor:
            labels_merged[i] = max_label
    return labels_merged


def order_labels(labels: list[int]) -> list[int]:
    labels_ordered = deepcopy(labels)
    num = 0
    d = defaultdict(int)
    for i, l in enumerate(labels_ordered):
        if not d[l]:
            num += 1
            d[l] = num
        labels_ordered[i] = d[l]
    return labels_ordered

