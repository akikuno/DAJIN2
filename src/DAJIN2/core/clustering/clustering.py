from __future__ import annotations

from collections import Counter
from itertools import chain
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix, spmatrix
from sklearn import metrics
from sklearn.cluster import BisectingKMeans

from DAJIN2.core.clustering.constrained_kmeans import ConstrainedKMeans
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
    positive_silhouette_previous = 1

    for i in range(2, coverage_sample):
        np.random.seed(seed=1)
        labels_all = (
            ConstrainedKMeans(n_clusters=i, min_cluster_size=min_cluster_size, random_state=1).fit_predict(X).tolist()
        )

        labels_sample = labels_all[:coverage_sample]
        labels_control = labels_all[coverage_sample:]

        min_cluster_size_ = min(Counter(labels_sample).values())

        num_labels_control = count_number_of_clusters(labels_control, coverage_control)

        silhouette_scores = metrics.silhouette_samples(X, labels_all, metric="euclidean")
        positive_silhouette_current = sum(silhouette_scores > 0.25)
        silhouette_ratio = positive_silhouette_current / positive_silhouette_previous

        # print(
        #     i,
        #     min_cluster_size,
        #     Counter(labels_control),
        #     Counter(labels_sample),
        #     silhouette_ratio,
        #     positive_silhouette_current,
        #     positive_silhouette_previous,
        # )  # ! DEBUG

        if min_cluster_size_ < min_cluster_size - 1 or num_labels_control >= 2 or silhouette_ratio < 0.95:
            break

        labels_previous = labels_sample
        positive_silhouette_previous = max(positive_silhouette_previous, positive_silhouette_current)
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
    X = np.array(list(chain(scores_sample, scores_control_subset)))

    coverage_sample = io.count_newlines(path_score_sample)
    coverage_control = len(scores_control_subset)
    labels = optimize_labels(X, coverage_sample, coverage_control, min_cluster_size)

    """Re-allocate clusters with strand bias to clusters without strand bias"""
    if strand_bias_in_control is False:
        labels = remove_biased_clusters(path_sample, path_score_sample, labels)

    return labels
