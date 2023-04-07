from __future__ import annotations

import warnings
from collections import Counter

import numpy as np
from sklearn.decomposition import PCA
from sklearn.exceptions import ConvergenceWarning
from sklearn.mixture import GaussianMixture

from DAJIN2.core.clustering.merge_clusters import merge_clusters

warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=ConvergenceWarning)

###############################################################################
# Dimension reduction
###############################################################################


def reduce_dimension(scores_sample: list[list], scores_control_subset: list[list]) -> np.array:
    scores = scores_sample + scores_control_subset
    n_components = min(20, len(scores[0]))
    pca = PCA(n_components=n_components).fit(scores)
    return pca.transform(scores)


def optimize_labels(X: np.array, scores_sample: list[list], scores_control_subset: list[list]) -> list[int]:
    scores = scores_sample + scores_control_subset
    n_components = min(20, len(scores))
    labels_results = [1] * len(scores_sample)
    for i in range(1, n_components):
        np.random.seed(seed=1)
        labels = GaussianMixture(n_components=i, random_state=1).fit_predict(X)
        labels = labels.tolist()
        labels_sample = labels[: len(scores_sample)]
        labels_control = labels[len(scores_sample) :]
        labels_merged = merge_clusters(labels_control, labels_sample)
        # Reads < 1% in the control are considered clustering errors and are not counted
        count_control = Counter(labels_control)
        num_labels_control = sum(1 for reads in count_control.values() if reads / sum(count_control.values()) * 100 > 1)
        if num_labels_control > 1:
            return labels_results
        labels_results = labels_merged
    return labels_results


###############################################################################
# main
###############################################################################


def return_labels(scores_sample: list[list], scores_control: list[list]) -> list[int]:
    np.random.seed(seed=1)
    X_control = reduce_dimension([], scores_control)
    labels = GaussianMixture(n_components=2, random_state=1).fit_predict(X_control)
    label_most = Counter(labels).most_common()[0][0]
    scores_control_subset = [s for l, s in zip(labels, scores_control) if l == label_most][:1000]
    X = reduce_dimension(scores_sample, scores_control_subset)
    labels = optimize_labels(X, scores_sample, scores_control_subset)
    return labels
