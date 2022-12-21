from __future__ import annotations
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.mixture import GaussianMixture
from collections import Counter

# from sklearn.cluster import MeanShift
from sklearn.decomposition import PCA

from src.DAJIN2.core.clustering.merge_clusters import merge_clusters
from src.DAJIN2.core.clustering.reorder_labels import reorder_labels

###############################################################################
# Dimention reduction
###############################################################################


def reduce_dimention(
    scores_sample: list[list[float]], scores_control_subset: list[list[float]], n_components: int = 20
) -> np.array:
    scores = scores_sample + scores_control_subset
    n_components = min(n_components, len(scores))
    scaler = StandardScaler()
    scores_scaler = scaler.fit_transform(scores)
    pca = PCA(n_components=n_components).fit(scores_scaler)
    variance = pca.explained_variance_
    return pca.transform(scores_scaler) * variance


def optimize_labels(
    X: np.array, scores_sample: list[list[float]], scores_control_subset: list[list[float]], n_components: int = 20
) -> list[int]:
    scores = scores_sample + scores_control_subset
    point_coodinates = []
    n_components = min(n_components, len(scores))
    for i in range(1, n_components):
        labels = GaussianMixture(n_components=i, random_state=0).fit_predict(X)
        labels = labels.tolist()
        labels_control = labels[len(scores_sample) :]
        labels_sample = labels[: len(scores_sample)]
        labels_merged = merge_clusters(labels_control, labels_sample)
        labels_reorder = reorder_labels(labels_merged)
        x = len(Counter(labels_control))
        y = len(Counter(labels_reorder))
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


def edist(x1, y1, x2, y2):
    return ((x1 - x2) ** 2 + (y1 - y2) ** 2) ** 0.5


def return_labels(
    scores_sample: list[list[float]], scores_control_subset: list[list[float]], n_components: int = 20
) -> list[int]:
    X = reduce_dimention(scores_sample, scores_control_subset)
    labels = optimize_labels(X, scores_sample, scores_control_subset)
    return labels
