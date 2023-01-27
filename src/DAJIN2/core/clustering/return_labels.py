from __future__ import annotations
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.mixture import GaussianMixture
from collections import Counter

from src.DAJIN2.core.clustering.merge_clusters import merge_clusters
from src.DAJIN2.core.clustering.reorder_labels import reorder_labels

###############################################################################
# Dimention reduction
###############################################################################


def reduce_dimention(scores_sample: list[list], scores_control_subset: list[list]) -> np.array:
    scores = scores_sample + scores_control_subset
    n_components = min(20, len(scores[0]))
    # scaler = StandardScaler()
    # scores = scaler.fit_transform(scores)
    pca = PCA(n_components=n_components).fit(scores)
    # variance = pca.explained_variance_
    return pca.transform(scores)  # * variance


def edist(x1, y1, x2, y2):
    return ((x1 - x2) ** 2 + (y1 - y2) ** 2) ** 0.5


def optimize_labels(X: np.array, scores_sample: list[list], scores_control_subset: list[list]) -> list[int]:
    scores = scores_sample + scores_control_subset
    point_coodinates = []
    n_components = min(20, len(scores))
    for i in range(1, n_components):
        np.random.seed(seed=1)
        labels = GaussianMixture(n_components=i, random_state=1).fit_predict(X)
        labels = labels.tolist()
        labels_sample = labels[: len(scores_sample)]
        labels_control = labels[len(scores_sample) :]
        labels_merged = merge_clusters(labels_control, labels_sample)
        labels_reorder = reorder_labels(labels_merged)
        x = len(Counter(labels_control))
        y = len(Counter(labels_reorder))
        print(i, Counter(labels_control), Counter(labels_reorder))  # ! -------------------------------------------
        point_coodinates.append([i, x, y, iter(labels_reorder)])
    idx = point_coodinates[0][0]
    x = max(c[1] for c in point_coodinates)
    y = max(c[2] for c in point_coodinates)
    dist_coodinates = [[i, edist(x1, y1, x, 0) / edist(x1, y1, 0, y)] for i, x1, y1, _ in point_coodinates]
    best_idx = sorted(dist_coodinates, key=lambda x: -x[1])[0][0]
    labels = list(point_coodinates[best_idx - idx][-1])
    return labels


###############################################################################
# main
###############################################################################


def return_labels(scores_sample: list[list], scores_control: list[list]) -> list[int]:
    np.random.seed(seed=1)
    X_control = reduce_dimention([], scores_control)
    labels = GaussianMixture(n_components=2, random_state=1).fit_predict(X_control)
    label_most = Counter(labels).most_common()[0][0]
    scores_control_subset = [s for l, s in zip(labels, scores_control) if l == label_most][:1000]
    X = reduce_dimention(scores_sample, scores_control_subset)
    labels = optimize_labels(X, scores_sample, scores_control_subset)
    return labels

