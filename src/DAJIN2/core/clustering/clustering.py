from __future__ import annotations

from DAJIN2.utils import io, config

config.set_warnings()

import numpy as np

from pathlib import Path
from itertools import chain
from typing import Generator
from collections import Counter
from collections import defaultdict

from sklearn import metrics
from sklearn.decomposition import PCA
from sklearn.mixture import GaussianMixture
from sklearn.tree import DecisionTreeClassifier

from DAJIN2.core.clustering.label_merger import merge_labels


###############################################################################
# Dimension reduction
###############################################################################


def reduce_dimension(scores_sample: Generator[list], scores_control: Generator[list]) -> np.array:
    scores = list(chain(scores_sample, scores_control))
    model = PCA(n_components=min(20, len(scores)))
    # return model.fit_transform(scores) * np.nan_to_num(model.explained_variance_ratio_)
    return model.fit_transform(scores)


###############################################################################
# Clustering
###############################################################################


def optimize_labels(X: np.array, coverage_sample, coverage_control) -> list[int]:
    n_components = min(20, coverage_sample + coverage_control)
    labels_prev = list(range(coverage_sample))
    for i in range(1, n_components):
        np.random.seed(seed=1)
        labels_all = GaussianMixture(n_components=i, random_state=1).fit_predict(X).tolist()
        labels_sample = labels_all[:coverage_sample]
        labels_control = labels_all[coverage_sample:]
        labels_merged = merge_labels(labels_control, labels_sample)
        # print(i, Counter(labels_sample), Counter(labels_control), Counter(labels_merged))  # ! DEBUG
        # Reads < 1% in the control are considered clustering errors and are not counted
        count_control = Counter(labels_control)
        num_labels_control = sum(1 for reads in count_control.values() if reads / coverage_control > 0.01)
        mutual_info = metrics.adjusted_mutual_info_score(labels_prev, labels_merged)
        # Report the number of clusters in SAMPLE when the number of clusters in CONTROL is split into more than one.
        if num_labels_control > 1 or 0.95 < mutual_info < 1.0:
            return labels_merged
        else:
            labels_results = labels_merged
        labels_prev = labels_merged
    return labels_results


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
        strand_biases[label] = not (0.25 < strand_ratio < 0.75)

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


def get_labels_removed_strand_bias(path_sample, path_score_sample, labels) -> list[int]:
    strand_biases = _get_strand_biases_on_each_label(labels, path_sample)
    # Until there is at least one True and one False or
    # 1000 iterations (1000 is a suitable number to exit an infinite loop just in case)
    i = 0
    while len(Counter(strand_biases.values())) > 1 and i < 1000:
        labels = _correct_clusters_with_strand_bias(path_score_sample, labels, strand_biases)
        strand_biases = _get_strand_biases_on_each_label(labels, path_sample)
        i += 1
    return labels
