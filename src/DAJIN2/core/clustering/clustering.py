from __future__ import annotations

from DAJIN2.utils import io, config

config.set_warnings()

# prevent BLAS from using all cores
config.set_single_threaded_blas()

from pathlib import Path
from itertools import chain

# from typing import Generator
from collections import Counter

import numpy as np
from sklearn import metrics
from sklearn.cluster import KMeans
from scipy.sparse import csr_matrix, spmatrix

# from sklearn.decomposition import PCA
# from sklearn.mixture import GaussianMixture

from DAJIN2.core.clustering.label_merger import merge_labels
from DAJIN2.core.clustering.score_handler import subset_scores
from DAJIN2.core.clustering.strand_bias_handler import remove_biased_clusters


###############################################################################
# Dimension reduction
###############################################################################


# def reduce_dimension(scores_sample: Generator[list], scores_control: Generator[list]) -> np.array:
#     scores = list(chain(scores_sample, scores_control))
#     model = PCA(n_components=min(20, len(scores)))
#     # return model.fit_transform(scores) * np.nan_to_num(model.explained_variance_ratio_)
#     return model.fit_transform(scores)


###############################################################################
# Clustering
###############################################################################


def optimize_labels(X: spmatrix, coverage_sample, coverage_control) -> list[int]:
    labels_results = list(range(coverage_sample))
    for i in range(1, coverage_sample):
        np.random.seed(seed=1)
        labels_all = KMeans(n_clusters=i, random_state=1, n_init="auto").fit_predict(X).tolist()
        labels_sample = labels_all[:coverage_sample]
        labels_control = labels_all[coverage_sample:]
        labels_merged = merge_labels(labels_control, labels_sample)
        # print(i, Counter(labels_sample), Counter(labels_control), Counter(labels_merged))  # ! DEBUG
        # Reads < 1% in the control are considered clustering errors and are not counted
        count_control = Counter(labels_control)
        num_labels_control = sum(1 for reads in count_control.values() if reads / coverage_control > 0.01)
        mutual_info = metrics.adjusted_rand_score(labels_results, labels_merged)
        # Report the number of clusters in SAMPLE when the number of clusters in CONTROL is split into more than one.
        if num_labels_control >= 2 or 0.95 <= mutual_info <= 1.0:
            return labels_results
        else:
            labels_results = labels_merged
    return labels_results


# def optimize_labels(X: np.array, coverage_sample, coverage_control) -> list[int]:
#     labels_prev = list(range(coverage_sample))
#     for i in range(1, coverage_sample):
#         np.random.seed(seed=1)
#         labels_all = GaussianMixture(n_components=i, random_state=1).fit_predict(X).tolist()
#         labels_sample = labels_all[:coverage_sample]
#         labels_control = labels_all[coverage_sample:]
#         labels_merged = merge_labels(labels_control, labels_sample)
#         # print(i, Counter(labels_sample), Counter(labels_control), Counter(labels_merged))  # ! DEBUG
#         # Reads < 1% in the control are considered clustering errors and are not counted
#         count_control = Counter(labels_control)
#         num_labels_control = sum(1 for reads in count_control.values() if reads / coverage_control > 0.01)
#         mutual_info = metrics.adjusted_mutual_info_score(labels_prev, labels_merged)
#         # Report the number of clusters in SAMPLE when the number of clusters in CONTROL is split into more than one.
#         if num_labels_control > 1 or 0.95 < mutual_info < 1.0:
#             return labels_prev
#         else:
#             labels_results = labels_merged
#         labels_prev = labels_merged
#     return labels_results


def get_label_most_common(labels: list[int]) -> int:
    return Counter(labels).most_common()[0][0]


def return_labels(path_score_sample: Path, path_score_control: Path, path_sample: Path, strand_bias: bool) -> list[int]:
    np.random.seed(seed=1)
    score_control = list(io.read_jsonl(path_score_control))
    X_control = csr_matrix(score_control)
    # subset to 1000 reads of controls in the most common cluster to remove outliers and reduce computation time
    labels_control = KMeans(n_clusters=2, random_state=1, n_init="auto").fit_predict(X_control)
    label_most_common = get_label_most_common(labels_control)
    scores_control_subset = subset_scores(labels_control, io.read_jsonl(path_score_control), label_most_common, 1000)
    scores_sample = list(io.read_jsonl(path_score_sample))
    X = csr_matrix(list(chain(scores_sample, scores_control_subset)))
    coverage_sample = io.count_newlines(path_score_sample)
    coverage_control = len(scores_control_subset)
    labels = optimize_labels(X, coverage_sample, coverage_control)
    # correct clusters with strand bias
    if strand_bias is False:
        labels = remove_biased_clusters(path_sample, path_score_sample, labels)
    return labels


# def return_labels(
#     path_score_sample: Path | str, path_score_control: Path | str, path_sample: Path | str, strand_bias: bool
# ) -> list[int]:
#     np.random.seed(seed=1)
#     X_control = reduce_dimension([], io.read_jsonl(path_score_control))
#     # subset to 1000 reads of controls in the most common cluster to remove outliers and reduce computation time
#     labels_control = GaussianMixture(n_components=2, random_state=1).fit_predict(X_control)
#     label_most_common = get_label_most_common(labels_control)
#     scores_control_subset = subset_scores(labels_control, io.read_jsonl(path_score_control), label_most_common, 1000)
#     X = reduce_dimension(io.read_jsonl(path_score_sample), scores_control_subset)
#     coverage_sample = io.count_newlines(path_score_sample)
#     coverage_control = len(scores_control_subset)
#     labels = optimize_labels(X, coverage_sample, coverage_control)
#     # correct clusters with strand bias
#     if strand_bias is False:
#         labels = remove_biased_clusters(path_sample, path_score_sample, labels)
#     return labels
