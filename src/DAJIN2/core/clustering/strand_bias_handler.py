from __future__ import annotations

from pathlib import Path
from collections import defaultdict, Counter

from DAJIN2.utils import io, config

# prevent BLAS from using all cores
config.set_single_threaded_blas()

from sklearn.tree import DecisionTreeClassifier

# Constants
STRAND_BIAS_LOWER_LIMIT = 0.1
STRAND_BIAS_UPPER_LIMIT = 0.9


def is_strand_bias(path_control: Path) -> bool:
    count_strand = defaultdict(int)
    for m in io.read_jsonl(path_control):
        count_strand[m["STRAND"]] += 1

    total = count_strand["+"] + count_strand["-"]
    percentage_plus = count_strand["+"] / total if total else 0

    return not (STRAND_BIAS_LOWER_LIMIT < percentage_plus < STRAND_BIAS_UPPER_LIMIT)


###############################################################################
# Handle Strand bias
# # Clusters of reads with mutations with strand bias are merged into similar clusters without strand bias
###############################################################################


def _count_strand(labels: list[int], samples: list[dict[str, str]]) -> tuple[defaultdict, defaultdict]:
    """Count the occurrences of each strand type by label."""
    count_strand_by_labels = defaultdict(int)
    total_count_by_labels = defaultdict(int)

    for label, sample in zip(labels, samples):
        total_count_by_labels[label] += 1
        if sample["STRAND"] == "+":
            count_strand_by_labels[label] += 1

    return count_strand_by_labels, total_count_by_labels


def _calculate_strand_biases(
    count_strand_by_labels: defaultdict, total_count_by_labels: defaultdict
) -> dict[int, bool]:
    """Calculate strand biases based on strand counts."""
    strand_biases = {}
    for label, total in total_count_by_labels.items():
        strand_count = count_strand_by_labels[label]
        strand_ratio = strand_count / total
        strand_biases[label] = not (STRAND_BIAS_LOWER_LIMIT < strand_ratio < STRAND_BIAS_UPPER_LIMIT)

    return strand_biases


def _get_strand_biases_on_each_label(labels: list[int], path_sample: Path | str) -> dict[int, bool]:
    """Get strand biases for given labels and samples.
    Args:
        labels: A list of integer labels.
        path_sample: The path to the sample file.
    Returns:
        A dictionary containing strand biases by label.
    """
    samples = io.read_jsonl(path_sample)
    count_strand_by_labels, total_count_by_labels = _count_strand(labels, samples)
    return _calculate_strand_biases(count_strand_by_labels, total_count_by_labels)


def _prepare_training_testing_sets(labels, scores, strand_biases) -> tuple[list, list, list]:
    x_train, y_train, x_test = [], [], []
    for label, score in zip(labels, scores):
        if strand_biases[label]:
            x_test.append(score)
        else:
            x_train.append(score)
            y_train.append(label)
    return x_train, y_train, x_test


def _train_decision_tree(x_train, y_train) -> DecisionTreeClassifier:
    dtree = DecisionTreeClassifier(random_state=1)
    dtree.fit(x_train, y_train)
    return dtree


def _allocate_labels(labels, strand_biases, dtree, x_test) -> list[int]:
    label_predictions = dtree.predict(x_test)
    label_predict_iter = iter(label_predictions)
    for i, label in enumerate(labels):
        if strand_biases[label]:
            labels[i] = next(label_predict_iter)
    return labels


def _correct_clusters_with_strand_bias(path_score_sample, labels, strand_biases) -> list[int]:
    scores = io.read_jsonl(path_score_sample)
    x_train, y_train, x_test = _prepare_training_testing_sets(labels, scores, strand_biases)
    dtree = _train_decision_tree(x_train, y_train)
    return _allocate_labels(labels, strand_biases, dtree, x_test)


def remove_biased_clusters(path_sample, path_score_sample, labels) -> list[int]:
    strand_biases = _get_strand_biases_on_each_label(labels, path_sample)
    # Until there is at least one True and one False or
    # 1000 iterations (1000 is a suitable number to exit an infinite loop just in case)
    i = 0
    labels_corrected = labels
    while len(Counter(strand_biases.values())) > 1 and i < 1000:
        labels_corrected = _correct_clusters_with_strand_bias(path_score_sample, labels_corrected, strand_biases)
        strand_biases = _get_strand_biases_on_each_label(labels_corrected, path_sample)
        i += 1
    return labels_corrected
