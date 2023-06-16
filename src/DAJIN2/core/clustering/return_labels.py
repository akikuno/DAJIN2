from __future__ import annotations

import warnings
from collections import Counter
from typing import Generator
from pathlib import Path
import json
from itertools import chain

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


def reduce_dimension(scores_sample: Generator[list], scores_control: Generator[list]) -> np.array:
    scores = list(chain(scores_sample, scores_control))
    model = PCA(n_components=min(20, len(scores)))
    return model.fit_transform(scores)


def optimize_labels(X: np.array, coverage_sample, coverage_control) -> list[int]:
    n_components = min(20, coverage_sample + coverage_control)
    for i in range(1, n_components):
        np.random.seed(seed=1)
        labels_all = GaussianMixture(n_components=i, random_state=1).fit_predict(X).tolist()
        labels_sample = labels_all[:coverage_sample]
        labels_control = labels_all[coverage_sample:]
        labels_merged = merge_clusters(labels_control, labels_sample)
        # print(i, Counter(labels_sample), Counter(labels_control), Counter(labels_merged))  # ! DEBUG
        # Reads < 1% in the control are considered clustering errors and are not counted
        count_control = Counter(labels_control)
        num_labels_control = sum(1 for reads in count_control.values() if reads / coverage_control > 0.01)
        # Report the number of clusters in SAMPLE when the number of clusters in CONTROL is split into more than one.
        if num_labels_control > 1:
            return labels_merged
        else:
            labels_results = labels_merged
    return labels_results


def read_json(filepath: Path | str) -> Generator[dict[str, str]]:
    with open(filepath, "r") as f:
        for line in f:
            yield json.loads(line)


def count_newlines(filepath: Path | str) -> int:
    def _make_gen(reader):
        while True:
            b = reader(2**16)
            if not b:
                break
            yield b

    with open(filepath, "rb") as f:
        count = sum(buf.count(b"\n") for buf in _make_gen(f.raw.read))
    return count


###############################################################################
# main
###############################################################################


def get_label_most(labels: list[int]) -> int:
    return Counter(labels).most_common()[0][0]


def subset_scores(labels: list[int], scores: list[int], label_most: int, length: int = 1000) -> list[int]:
    subset = [score for label, score in zip(labels, scores) if label == label_most]
    return subset[:length]


def return_labels(path_score_sample: Path | str, path_score_control: Path | str) -> list[int]:
    np.random.seed(seed=1)
    X_control = reduce_dimension([], read_json(path_score_control))
    # subset to 1000 reads of controls in the most common cluster to remove outliers and reduce computation time
    labels_control = GaussianMixture(n_components=2, random_state=1).fit_predict(X_control)
    label_most = get_label_most(labels_control)
    scores_control_subset = subset_scores(labels_control, read_json(path_score_control), label_most, 1000)
    X = reduce_dimension(read_json(path_score_sample), scores_control_subset)
    coverage_sample = count_newlines(path_score_sample)
    coverage_control = len(scores_control_subset)
    labels = optimize_labels(X, coverage_sample, coverage_control)
    return labels
