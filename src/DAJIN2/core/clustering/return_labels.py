import numpy as np
import scipy.stats as st
from sklearn.cluster import Birch
from scipy.sparse import csr_matrix
from collections import defaultdict


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


def return_labels(scores):
    # flatten from 3d to 2d
    xarr = np.array(scores)
    xflatten = np.reshape(xarr, (len(xarr), -1))
    xflatten = csr_matrix(xflatten)
    n_clusters = 1
    while True or n_clusters < 1000:
        prev_labels = predict_labels(xflatten, n_clusters=n_clusters)
        prev_table = make_table(prev_labels)
        current_labels = predict_labels(xflatten, n_clusters=n_clusters + 1)
        current_table = make_table(current_labels)
        if n_clusters == 1:
            prev_table += [0]
        else:
            current_table = current_table[:-1]
        pval = chistatistic(prev_table, current_table, threshold=0.05)
        if pval > 0.05:
            labels = prev_labels
            break
        n_clusters += 1
    return labels

