from __future__ import annotations

from collections import Counter


def calculate_label_percentages(labels: list[int]) -> dict[int, float]:
    """Calculate the percentage of each label in the list of labels."""
    total_labels = len(labels)
    label_counts = Counter(labels)
    return {label: (count / total_labels * 100) for label, count in label_counts.items()}


def merge_mixed_cluster(labels_control: list[int], labels_sample: list[int], threshold: float = 0.5) -> list[int]:
    """Merge labels in sample if they appear more than 'threshold' percentage in control."""
    labels_merged = labels_sample.copy()
    label_percentages_control = calculate_label_percentages(labels_control)
    mixed_labels = {label for label, percent in label_percentages_control.items() if percent > threshold}

    new_label = max(labels_merged) + 1
    for i, label in enumerate(labels_sample):
        if label in mixed_labels:
            labels_merged[i] = new_label

    return labels_merged


def merge_minor_cluster(labels_sample: list[int], threshold: float = 0.5) -> list[int]:
    """Merge labels in sample if they appear less than 'threshold' percentage."""
    labels_merged = labels_sample.copy()
    label_percentages_sample = calculate_label_percentages(labels_sample)
    minor_labels = {label for label, percent in label_percentages_sample.items() if percent < threshold}

    most_common_label = Counter(labels_sample).most_common(1)[0][0]
    for i, label in enumerate(labels_sample):
        if label in minor_labels:
            labels_merged[i] = most_common_label

    return labels_merged


###############################################################################
# main
###############################################################################


def merge_labels(labels_control: list[int], labels_sample: list[int]) -> list[int]:
    labels_merged = merge_mixed_cluster(labels_control, labels_sample)
    labels_merged = merge_minor_cluster(labels_merged)
    return labels_merged
