from __future__ import annotations
from sklearn.cluster import MeanShift

###############################################################################
# main
###############################################################################


def return_labels(scores: list[list[float]], THREADS: int = 1) -> list[int]:
    labels = MeanShift(n_jobs=THREADS).fit(scores).labels_
    return labels.tolist()
