from __future__ import annotations

from collections import defaultdict
from pathlib import Path

from scipy.stats import fisher_exact
from sklearn.tree import DecisionTreeClassifier

from DAJIN2.utils import io

"""
Nanopore sequencing results often results in strand specific mutations even though the mutation is not strand specific, thus they are considered as sequencing errors and should be removed.

This module provides functions to determine whether each allele obtained after clustering is formed due to sequencing errors caused by strand bias.

Re-allocates reads belonging to clusters with strand bias to clusters without strand bias.
"""


def count_strand_distribution(path_midsv: Path) -> dict[str, int]:
    """
    Counts the distribution of '+' and '-' strands in the given control file.
    """
    strand_count = defaultdict(int)
    for sample in io.read_jsonl(path_midsv):
        strand_count[sample["STRAND"]] += 1

    # Ensure '+' and '-' are present for all labels
    strand_count.setdefault("+", 0)
    strand_count.setdefault("-", 0)

    return dict(strand_count)


###############################################################################
# Handle Strand bias
# # Clusters of reads with mutations with strand bias are merged into similar clusters without strand bias
###############################################################################


def count_strand_by_label(samples: dict[str, str], labels: list[int]) -> dict[int, dict[str, int]]:
    """Count the occurrences of each strand type ('+' and '-') for each label."""
    strand_counts_by_label = defaultdict(lambda: defaultdict(int))

    for label, sample in zip(labels, samples):
        strand_counts_by_label[label][sample["STRAND"]] += 1

    # Ensure '+' and '-' are present for all labels
    for counts in strand_counts_by_label.values():
        counts.setdefault("+", 0)
        counts.setdefault("-", 0)

    return {label: dict(counts) for label, counts in strand_counts_by_label.items()}


def determine_strand_biases(
    strand_counts_by_labels: dict[int, dict[str, int]], strand_sample: dict[str, int]
) -> dict[int, bool]:
    """Determine strand biases based on positive strand counts."""
    strand_biases = {}
    for label, strand_counts in strand_counts_by_labels.items():
        # numerical bias detection
        plus_ratio = strand_counts["+"] / (strand_counts["+"] + strand_counts["-"])
        is_biased_ratio = False if 0.2 < plus_ratio < 0.8 else True
        # statistical bias detection
        table = [[strand_counts["+"], strand_counts["-"]], [strand_sample["+"], strand_sample["-"]]]
        res = fisher_exact(table, alternative="two-sided")

        strand_biases[label] = is_biased_ratio and res.pvalue < 0.01

    return strand_biases


def annotate_strand_bias_by_labels(path_sample: Path, labels: list[int], strand_sample: dict[str, int]) -> bool:
    """Determine whether there is strand bias in the samples based on the provided labels."""
    samples = io.read_jsonl(path_sample)
    strand_counts_by_labels = count_strand_by_label(samples, labels)
    return determine_strand_biases(strand_counts_by_labels, strand_sample)


def prepare_training_testing_sets(labels, scores, strand_biases) -> tuple[list, list, list]:
    """Prepare training and testing datasets based on strand biases."""
    train_data, train_labels, test_data = [], [], []
    for label, score in zip(labels, scores):
        if strand_biases[label]:
            test_data.append(score)
        else:
            train_data.append(score)
            train_labels.append(label)
    return train_data, train_labels, test_data


def train_decision_tree(train_data, train_labels) -> DecisionTreeClassifier:
    """Train a decision tree classifier using the provided features and labels."""
    dtree = DecisionTreeClassifier(random_state=1)
    dtree.fit(train_data, train_labels)
    return dtree


def allocate_labels(labels: list[int], strand_biases: dict[str, bool], dtree, test_data) -> list[int]:
    """Re-allocates reads belonging to clusters with strand bias to clusters without strand bias."""
    label_predictions = iter(dtree.predict(test_data))
    for i, label in enumerate(labels):
        if strand_biases[label]:
            labels[i] = next(label_predictions).tolist()
    return labels


def remove_biased_clusters(
    path_sample: Path, path_score_sample: Path, labels: list[int], strand_sample: dict[str, int]
) -> list[int]:
    """Remove clusters with strand bias by re-labeling based on decision tree predictions.
    Continue until at least one of the samples exhibits strand bias (i.e., do not calculate if all samples exhibit strand bias or, if none of the samples exhibit strand bias) and
    1000 iterations are reached, which serves as a safeguard to prevent infinite loops.
    """
    strand_biases = annotate_strand_bias_by_labels(path_sample, labels, strand_sample)

    iteration_count = 0
    labels_corrected = labels
    while len(set(strand_biases.values())) > 1 and iteration_count < 1000:
        # Re-allocation of labels of biased clusters to unbiased clusters
        scores = io.read_jsonl(path_score_sample)
        train_data, train_labels, test_data = prepare_training_testing_sets(labels, scores, strand_biases)
        dtree = train_decision_tree(train_data, train_labels)
        labels_corrected = allocate_labels(labels, strand_biases, dtree, test_data)

        # Re-calculate strand biases based on the corrected labels
        strand_biases = annotate_strand_bias_by_labels(path_sample, labels_corrected, strand_sample)

        iteration_count += 1

    return labels_corrected
