from __future__ import annotations

from collections import Counter
from itertools import chain
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix, spmatrix
from sklearn import metrics
from sklearn.cluster import BisectingKMeans

from DAJIN2.core.clustering.label_merger import merge_labels
from DAJIN2.core.clustering.score_handler import subset_scores
from DAJIN2.core.clustering.strand_bias_handler import remove_biased_clusters
from DAJIN2.utils import config, io

config.set_warnings_ignore()


###############################################################################
# Clustering
###############################################################################


def count_number_of_clusters(labels_control: list[int], coverage_control: int) -> int:
    """If there is less than 1% lead within a cluster, they are considered clustering errors and are not counted"""
    return sum(1 for control_reads in Counter(labels_control).values() if control_reads / coverage_control > 0.01)


def optimize_labels(X: spmatrix, coverage_sample: int, coverage_control: int, min_cluster_size: int) -> list[int]:
    labels_previous = [1 for _ in range(coverage_sample)]

    for i in range(2, coverage_sample):
        np.random.seed(seed=1)
        labels_all = BisectingKMeans(n_clusters=i, random_state=1).fit_predict(X).tolist()

        labels_sample = labels_all[:coverage_sample]
        labels_control = labels_all[coverage_sample:]
        labels_current = merge_labels(labels_control, labels_sample, labels_previous)
        # print(i, Counter(labels_sample), Counter(labels_control), Counter(labels_current))  # ! DEBUG

        num_labels_control = count_number_of_clusters(labels_control, coverage_control)
        rand_index = metrics.adjusted_rand_score(labels_previous, labels_current)

        """
        Return the number of clusters when:
            - the number of clusters in control is split into more than one.
            - the mutual information between the current and previous labels is high enough (= similar).
            - the minimum cluster size is reached.
        To reduce the allele number, previous labels are returned.
        """
        if num_labels_control >= 2 or rand_index >= 0.95 or min_cluster_size > min(Counter(labels_current).values()):
            return labels_previous

        labels_previous = labels_current

    return labels_previous


def get_label_most_common(labels: list[int]) -> int:
    return Counter(labels).most_common()[0][0]


def return_labels(
    path_score_sample: Path,
    path_score_control: Path,
    path_sample: Path,
    strand_bias_in_control: bool,
    min_cluster_size: int,
) -> list[int]:
    np.random.seed(seed=1)

    score_control = list(io.read_jsonl(path_score_control))
    X_control = csr_matrix(score_control)

    """Subset to 1000 reads of controls in the most common cluster to remove outliers and reduce computation time"""
    labels_control = BisectingKMeans(n_clusters=2, random_state=1).fit_predict(X_control)
    label_most_common = get_label_most_common(labels_control)
    scores_control_subset = subset_scores(labels_control, io.read_jsonl(path_score_control), label_most_common, 1000)

    scores_sample = list(io.read_jsonl(path_score_sample))
    X = csr_matrix(list(chain(scores_sample, scores_control_subset)))

    coverage_sample = io.count_newlines(path_score_sample)
    coverage_control = len(scores_control_subset)
    labels = optimize_labels(X, coverage_sample, coverage_control, min_cluster_size)

    """Re-allocate clusters with strand bias to clusters without strand bias"""
    if strand_bias_in_control is False:
        labels = remove_biased_clusters(path_sample, path_score_sample, labels)

    return labels
