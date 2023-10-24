from __future__ import annotations


def add_labels(classif_sample: list[dict], labels: list[int]) -> list[dict]:
    """Add 'LABEL' key to each sample dict, indicating each clster of allele."""
    clust_sample = classif_sample.copy()

    for clust, label in zip(clust_sample, labels):
        clust["LABEL"] = label

    return clust_sample


def count_labels(clust_sample: list[dict]) -> dict[int, int]:
    """Count the occurrences of each label in the cluster sample."""
    label_count = {}
    for sample in clust_sample:
        label = sample["LABEL"]
        label_count[label] = label_count.get(label, 0) + 1
    return label_count


def add_readnum(clust_sample: list[dict]) -> list[dict]:
    """Add 'READNUM' key to each sample dict, indicating the number of occurrences of each label."""
    clust_result = clust_sample.copy()
    label_count = count_labels(clust_result)

    for sample in clust_result:
        sample["READNUM"] = label_count[sample["LABEL"]]

    return clust_result


def add_percent(clust_sample: list[dict]) -> list[dict]:
    """Add 'PERCENT' key to each sample dict, indicating the percentage of occurrences of each label."""
    clust_result = clust_sample.copy()
    total_samples = len(clust_result)
    label_count = count_labels(clust_result)

    for sample in clust_result:
        label = sample["LABEL"]
        sample["PERCENT"] = round((label_count[label] / total_samples) * 100, 3)

    return clust_result
