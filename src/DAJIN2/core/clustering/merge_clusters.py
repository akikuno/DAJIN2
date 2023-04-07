from __future__ import annotations

from collections import Counter, defaultdict


def merge_mixed_cluster(labels_control: list[int], labels_sample: list[int]) -> list[int]:
    labels_merged = labels_sample.copy()
    coverage_control = len(labels_control)
    labels_percent_control = defaultdict(float)
    for i, label in enumerate(labels_control):
        labels_percent_control[label] += 1 / coverage_control * 100
    labels_mixed = {i for i, per in labels_percent_control.items() if per > 0.5}
    max_label = max(labels_merged) + 1
    for i, label in enumerate(labels_sample):
        if label in labels_mixed:
            labels_merged[i] = max_label
    return labels_merged


def merge_minor_cluster(labels_sample: list[int]) -> list[int]:
    total_num = sum(Counter(labels_sample).values())
    cnt = Counter(labels_sample)
    labels_minor = {label for label, num in cnt.items() if num / total_num * 100 < 0.5}
    comon_label = cnt.most_common()[0][0]
    labels_merged = labels_sample.copy()
    for i, label in enumerate(labels_sample):
        if label in labels_minor:
            labels_merged[i] = comon_label
    return labels_merged


###############################################################################
# main
###############################################################################


def merge_clusters(labels_control: list[int], labels_sample: list[int]) -> list[int]:
    labels_merged = merge_mixed_cluster(labels_control, labels_sample)
    labels_merged = merge_minor_cluster(labels_merged)
    return labels_merged
