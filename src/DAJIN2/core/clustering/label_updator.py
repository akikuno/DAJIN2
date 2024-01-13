from __future__ import annotations


def relabel_with_consective_order(labels: list[int], start: int = 1) -> list[int]:
    """
    Reorders the given list of labels starting from the given start value.

    Args:
    - labels: List of integer labels to be reordered
    - start: The starting value for reordering

    Returns:
    - List of reordered labels
    """
    reordered_labels = []
    next_value = start
    label_to_value = {}  # Mapping from original label to new value

    for label in labels:
        if label not in label_to_value:
            label_to_value[label] = next_value
            next_value += 1

        reordered_labels.append(label_to_value[label])

    return reordered_labels


def update_labels(clust_sample: list[dict]) -> list[dict]:
    """
    Allocate new labels according to the ranking by PERCENT
    """
    clust_result = clust_sample.copy()
    clust_result.sort(key=lambda x: (-x["PERCENT"], x["LABEL"]))
    new_label = 1
    prev_label = clust_result[0]["LABEL"]
    for cs in clust_result:
        if prev_label != cs["LABEL"]:
            new_label += 1
        prev_label = cs["LABEL"]
        cs["LABEL"] = new_label
    return clust_result
