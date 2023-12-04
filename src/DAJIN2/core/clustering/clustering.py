from __future__ import annotations

from DAJIN2.utils import io, config

config.set_warnings()

# prevent BLAS from using all cores
config.set_single_threaded_blas()

from pathlib import Path
from itertools import chain

from collections import Counter

import numpy as np
from sklearn import metrics
from sklearn.cluster import KMeans
from scipy.sparse import csr_matrix, spmatrix

from DAJIN2.core.clustering.label_merger import merge_labels
from DAJIN2.core.clustering.score_handler import subset_scores
from DAJIN2.core.clustering.strand_bias_handler import remove_biased_clusters


###############################################################################
# Clustering
###############################################################################


def optimize_labels(X: spmatrix, coverage_sample, coverage_control) -> list[int]:
    labels_previous = list(range(coverage_sample))
    for i in range(1, coverage_sample):
        np.random.seed(seed=1)
        labels_all = KMeans(n_clusters=i, random_state=1, n_init="auto").fit_predict(X).tolist()
        labels_sample = labels_all[:coverage_sample]
        labels_control = labels_all[coverage_sample:]
        labels_current = merge_labels(labels_control, labels_sample)
        # print(i, Counter(labels_sample), Counter(labels_control), Counter(labels_current))  # ! DEBUG
        # Reads < 1% in the control are considered clustering errors and are not counted
        count_control = Counter(labels_control)
        num_labels_control = sum(1 for reads in count_control.values() if reads / coverage_control > 0.01)
        mutual_info = metrics.adjusted_rand_score(labels_previous, labels_current)
        # Report the number of clusters in SAMPLE when the number of clusters in CONTROL is split into more than one.
        if num_labels_control >= 2:
            return labels_previous
        if 0.95 <= mutual_info <= 1.0:
            return labels_current
        labels_previous = labels_current
    return labels_previous


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
