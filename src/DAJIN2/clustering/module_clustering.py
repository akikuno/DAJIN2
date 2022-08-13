import numpy as np
import scipy.stats as st
from sklearn.cluster import Birch
from scipy.sparse import csr_matrix
from collections import defaultdict


def return_labels(xflatten, n_clusters=1):
    brc = Birch(n_clusters=n_clusters)
    brc.fit(xflatten)
    return brc.predict(xflatten)


def make_table(labels, n_clusters=1):
    count = defaultdict(int)
    for label in labels:
        count[label] += 1
    table = list(count.values())
    table.sort(reverse=True)
    if n_clusters == 1:
        table += [0]
    elif n_clusters == 2:
        pass
    else:
        table = table[:-1]
    return table


def chistatistic(s_table, c_table, threshold=0.05) -> float:
    chisq_value, _, ddof, _ = st.chi2_contingency([s_table, c_table])
    delta = sum(s_table + c_table) * threshold ** 2
    pval = 1 - st.ncx2.cdf(chisq_value, ddof, delta)
    return pval


def clustering(scores):
    # flatten from 3d to 2d
    xarr = np.array(scores)
    xflatten = np.reshape(xarr, (len(xarr), -1))
    xflatten = csr_matrix(xflatten)
    # Clustering
    prev_labels = [0 for _ in range(len(scores))]
    prev_labels = np.array(prev_labels)
    prev_table = make_table(prev_labels, n_clusters=1)
    n_clusters = 2
    while True:
        current_labels = return_labels(xflatten, n_clusters=n_clusters)
        current_table = make_table(current_labels, n_clusters=n_clusters)
        pval = chistatistic(prev_table, current_table, threshold=0.05)
        if pval > 0.05:
            labels = prev_labels
            break
        prev_table = current_table
        prev_labels = current_labels
        n_clusters += 1
    return labels

