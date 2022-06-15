import sys
import os
import numpy as np
import pandas as pd
import hdbscan
import joblib
from sklearn.decomposition import PCA

################################################################################
# Input
################################################################################

args = sys.argv
if args != [""]:
    df_score = pd.read_csv(args[1], header=None)
    threads = int(args[2])
else:
    df_score = pd.read_csv(".DAJIN_temp/clustering/tmp_score.csv", header=None)
    threads = os.cpu_count() - 1

tmp = (df_score.sum(axis=0) != 0).to_list()
df_score = df_score.loc[:, tmp]

################################################################################
# PCA as a preprocessing step
################################################################################

pca = PCA(n_components=5)
pc_score = pca.fit_transform(df_score)
prop_var = pca.explained_variance_ratio_
pc_score = pd.DataFrame(pc_score) * prop_var

################################################################################
# Clustering
################################################################################

# ===========================================================
# Extract cluster size with the largest cluster size
# and the most frequent cluster numbers
# ===========================================================


def hdb_cl_num(cl_size):
    cl = hdbscan.HDBSCAN(
        min_samples=1,
        min_cluster_size=cl_size,
        core_dist_n_jobs=threads,
        memory=joblib.Memory(location=".DAJIN_temp/clustering", verbose=0),
    )
    tmp = cl.fit_predict(pc_score)
    return len(np.unique(tmp))


def hdb_predict(data, cl_size):
    cl = hdbscan.HDBSCAN(
        min_samples=1,
        min_cluster_size=cl_size,
        core_dist_n_jobs=threads,
        memory=joblib.Memory(location=".DAJIN_temp/clustering", verbose=0),
    )
    return cl.fit_predict(data) + 1


cl_range = [int(n) for n in np.linspace(
    len(pc_score) * 0.2, len(pc_score) * 0.4, 10) + 2]

hdb = [hdb_cl_num(cl_num) for cl_num in cl_range]
# hdb = [i for i in hdb if i != 1]

if len(hdb) > 0:
    _uniqs, _counts = np.unique(hdb, return_counts=True)
    max_mode = max(_uniqs[_counts == np.amax(_counts)])
    cl_index = max([i for i, x in enumerate(hdb) if x == max_mode])
    cl_size = cl_range[cl_index]
else:
    cl_size = max(cl_range)


cl_hdb = hdb_predict(pc_score, cl_size) + 2

# ===========================================================
# Repeat clustering
# ===========================================================

for _cl in np.unique(cl_hdb):
    if sum(cl_hdb == _cl) <= 5:
        continue
    _df_score = df_score[cl_hdb == _cl]
    _pca = PCA(n_components=5)
    _pc_score = _pca.fit_transform(_df_score)
    _prop_var = _pca.explained_variance_ratio_
    _pc_score = pd.DataFrame(_pc_score) * _prop_var
    _cl_size = int(len(_pc_score) * 0.4 + 2)
    _cl_hdb = hdb_predict(_pc_score, _cl_size) + 2
    cl_hdb[cl_hdb == _cl] = _cl_hdb + (_cl * 100)

cl_hdb = pd.Series(cl_hdb).rank(method="dense").astype(int)

################################################################################
# Output
################################################################################

cl_hdb.to_csv(sys.stdout, index=False, header=False)
