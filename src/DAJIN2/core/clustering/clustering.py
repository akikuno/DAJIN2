from __future__ import annotations

from DAJIN2.utils import config

config.set_warnings()

from itertools import chain
from typing import Generator
from collections import Counter

import numpy as np
from sklearn import metrics
from sklearn.decomposition import PCA
from sklearn.mixture import GaussianMixture

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
    labels_prev = list(range(coverage_sample))
    for i in range(1, coverage_sample):
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
