from __future__ import annotations
import numpy as np
import scipy.stats as st
from sklearn.cluster import OPTICS
from sklearn.cluster import Birch
from collections import defaultdict

np.seterr(divide="ignore")  # to suppress "RuntimeWarning: divide by zero encountered in true_divide"
# from scipy.sparse import csr_matrix


def predict_labels(xflatten, n_clusters=1):
    brc = Birch(n_clusters=n_clusters)
    brc.fit(xflatten)
    return brc.predict(xflatten)


def make_table(labels):
    count = defaultdict(int)
    for label in labels:
        count[label] += 1
    table = list(count.values())
    table.sort(reverse=True)
    return table


def chistatistic(prev_table, current_table, threshold=0.05) -> float:
    chisq_value, _, ddof, _ = st.chi2_contingency([prev_table, current_table])
    delta = sum(prev_table + current_table) * threshold ** 2
    pval = 1 - st.ncx2.cdf(chisq_value, ddof, delta)
    return pval


def return_labels(scores: list[list[float]]) -> list[int]:
    scores_flatten = np.reshape(scores, (len(scores), -1))
    labels = OPTICS(min_samples=0.001).fit_predict(scores_flatten)
    labels[labels == -1] = np.max(labels) + 1
    return labels.tolist()


# def return_labels(scores) -> list[int]:
#     labels = [0] * len(scores)
#     scores_flatten = np.reshape(scores, (len(scores), -1))
#     scores_csr = csr_matrix(scores_flatten)
#     n_clusters = 2
#     while True or n_clusters < 1000:
#         current_labels = predict_labels(scores_csr, n_clusters=n_clusters)
#         if n_clusters > len(np.unique(current_labels)):  # For ConvergenceWarning in BIRCH
#             break
#         else:
#             silhouette = metrics.silhouette_score(scores_csr, current_labels, metric="euclidean")
#             if silhouette > 0.99:
#                 labels = current_labels
#         n_clusters += 1
#     return labels


# def return_labels(scores) -> list[int]:
#     # flatten from 3d to 2d
#     scores_flatten = np.reshape(scores, (len(scores), -1))
#     scores_csr = csr_matrix(scores_flatten)
#     n_clusters = 1
#     while True or n_clusters < 1000:
#         prev_labels = predict_labels(scores_csr, n_clusters=n_clusters)
#         prev_table = make_table(prev_labels)
#         current_labels = predict_labels(scores_csr, n_clusters=n_clusters + 1)
#         current_table = make_table(current_labels)
#         if len(prev_table) == len(current_table):  # For ConvergenceWarning in BIRCH
#             labels = prev_labels
#             break
#         elif n_clusters == 1:
#             prev_table += [0]
#         else:
#             current_table = current_table[:-1]
#         pval = chistatistic(prev_table, current_table, threshold=0.05)
#         if pval > 0.05:
#             labels = prev_labels
#             break
#         n_clusters += 1
#     return labels.tolist()

