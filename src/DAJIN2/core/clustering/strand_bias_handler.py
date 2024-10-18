from __future__ import annotations

from collections import defaultdict
from collections.abc import Iterator
from pathlib import Path

from sklearn.tree import DecisionTreeClassifier

from DAJIN2.utils import io

"""
Nanopore sequencing results often results in strand specific mutations even though the mutation is not strand specific, thus they are considered as sequencing errors and should be removed.

This module provides functions to determine whether each allele obtained after clustering is formed due to sequencing errors caused by strand bias.

Re-allocates reads belonging to clusters with strand bias to clusters without strand bias.
"""

# Constants
STRAND_BIAS_LOWER_LIMIT = 0.2
STRAND_BIAS_UPPER_LIMIT = 0.8


def is_strand_bias(path_control: Path) -> bool:
    """
    Determines whether there is a strand bias in sequencing data
    based on the distribution of '+' and '-' strands.
    """
    count_strand = defaultdict(int)
    for sample in io.read_jsonl(path_control):
        count_strand[sample["STRAND"]] += 1

    total = count_strand["+"] + count_strand["-"]
    percentage_plus = count_strand["+"] / total if total > 0 else 0

    return not (STRAND_BIAS_LOWER_LIMIT < percentage_plus < STRAND_BIAS_UPPER_LIMIT)


###############################################################################
# Handle Strand bias
# # Clusters of reads with mutations with strand bias are merged into similar clusters without strand bias
###############################################################################


def count_strand(labels: list[int], samples: Iterator[dict[str, str]]) -> tuple[dict[str, int], dict[str, int]]:
    """Count the occurrences of each strand type by label."""
    positive_strand_counts_by_labels = defaultdict(int)
    total_counts_by_labels = defaultdict(int)

    for label, sample in zip(labels, samples):
        total_counts_by_labels[label] += 1
        if sample["STRAND"] == "+":
            positive_strand_counts_by_labels[label] += 1

    return dict(positive_strand_counts_by_labels), dict(total_counts_by_labels)


def determine_strand_biases(
    positive_strand_counts_by_labels: dict[str, int], total_counts_by_labels: dict[str, int]
) -> dict[int, bool]:
    """Determine strand biases based on positive strand counts."""
    strand_biases = {}
    for label, total in total_counts_by_labels.items():
        positive_strand_count = positive_strand_counts_by_labels.get(label, 0)
        strand_ratio = positive_strand_count / total
        strand_biases[label] = not (STRAND_BIAS_LOWER_LIMIT < strand_ratio < STRAND_BIAS_UPPER_LIMIT)

    return strand_biases


def annotate_strand_bias_by_labels(path_sample: Path, labels: list[int]) -> bool:
    """Determine whether there is strand bias in the samples based on the provided labels."""
    samples = io.read_jsonl(path_sample)
    positive_strand_counts_by_labels, total_counts_by_labels = count_strand(labels, samples)
    return determine_strand_biases(positive_strand_counts_by_labels, total_counts_by_labels)


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
            labels[i] = next(label_predictions)
    return labels


def remove_biased_clusters(path_sample: Path, path_score_sample: Path, labels: list[int]) -> list[int]:
    """Remove clusters with strand bias by re-labeling based on decision tree predictions.
    Continue until at least one of the samples exhibits strand bias (i.e., do not calculate if all samples exhibit strand bias or, if none of the samples exhibit strand bias) and
    1000 iterations are reached, which serves as a safeguard to prevent infinite loops.
    """
    strand_biases = annotate_strand_bias_by_labels(path_sample, labels)

    iteration_count = 0
    labels_corrected = labels
    while len(set(strand_biases.values())) > 1 and iteration_count < 1000:
        # Re-allocation of labels of biased clusters to unbiased clusters
        scores = io.read_jsonl(path_score_sample)
        train_data, train_labels, test_data = prepare_training_testing_sets(labels, scores, strand_biases)
        dtree = train_decision_tree(train_data, train_labels)
        labels_corrected = allocate_labels(labels, strand_biases, dtree, test_data)

        # Re-calculate strand biases based on the corrected labels
        strand_biases = annotate_strand_bias_by_labels(path_sample, labels)

        iteration_count += 1

    return labels_corrected
