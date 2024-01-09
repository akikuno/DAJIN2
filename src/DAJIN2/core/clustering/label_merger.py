from __future__ import annotations

from collections import Counter
import numpy as np


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


def map_clusters_to_previous(labels_sample: list[int], labels_previous: list[int]) -> dict[int, int]:
    """
    Determine which cluster in labels_previous corresponds to each cluster in labels_sample.
    """
    correspondence = {}
    for current_label in np.unique(labels_sample):
        # Get the indices of data points belonging to a specific cluster in labels_sample
        indices = np.where(labels_sample == current_label)[0]
        # Identify the cluster in labels_previous that these data points belong to
        previous_label = np.bincount(np.array(labels_previous)[indices]).argmax()
        correspondence[current_label] = previous_label

    return correspondence


def merge_minor_cluster(
    labels_sample: list[int],
    labels_previous: list[int],
    threshold_percentage: float = 0.5,
    threshold_readnumber: int = 10,
) -> list[int]:
    """Merge labels in sample if they appear less than 'threshold' percentage."""

    # Find minor labels
    label_percentages = calculate_label_percentages(labels_sample)
    minor_labels_percentage = {label for label, percent in label_percentages.items() if percent < threshold_percentage}
    minor_labels_readnumber = {label for label, num in Counter(labels_sample).items() if num <= threshold_readnumber}
    minor_labels = minor_labels_percentage | minor_labels_readnumber

    correspondence = map_clusters_to_previous(labels_sample, labels_previous)

    labels_merged = [correspondence[label] if label in minor_labels else label for label in labels_sample]

    return labels_merged


###############################################################################
# main
###############################################################################


def merge_labels(labels_control: list[int], labels_sample: list[int], labels_previous: list[int]) -> list[int]:
    labels_merged = merge_mixed_cluster(labels_control, labels_sample)
    labels_merged = merge_minor_cluster(
        labels_merged, labels_previous, threshold_percentage=0.5, threshold_readnumber=10
    )
    return labels_merged
