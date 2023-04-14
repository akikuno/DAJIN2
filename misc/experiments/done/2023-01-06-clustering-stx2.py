import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import midsv
from DAJIN2.core.clustering.preprocess import replace_both_ends_n
from src.DAJIN2.core.clustering.make_score import make_score
from src.DAJIN2.core.clustering.annotate_score import annotate_score
from src.DAJIN2.core.clustering.reorder_labels import reorder_labels
from src.DAJIN2.core.clustering.return_labels import return_labels
from collections import *

from sklearn.preprocessing import StandardScaler
from sklearn.mixture import GaussianMixture
from collections import Counter

# from sklearn.cluster import MeanShift
from sklearn.decomposition import PCA

from src.DAJIN2.core.clustering.merge_clusters import merge_clusters
from src.DAJIN2.core.clustering.reorder_labels import reorder_labels


small = "836991e0ff07"
medium = "d2a68de457cc"
large = "a0d18c462d05"


allele = "deletion"
sv = True
midsv_sample = midsv.read_jsonl("DAJINResults/.tempdir/test-knockout/midsv/test_barcode25_splice_deletion.jsonl")
midsv_control = midsv.read_jsonl("DAJINResults/.tempdir/test-knockout/midsv/test_barcode30_splice_deletion.jsonl")


midsv_test = []
for samp in midsv_sample:
    if any(s in samp["QNAME"] for s in [small, medium, large]):
        midsv_test.append(samp)

len(midsv_test)

# cssplits_sample = [cs["CSSPLIT"].split(",") for cs in midsv_test]
cssplits_sample = [cs["CSSPLIT"].split(",") for cs in midsv_sample]
cssplits_control = [cs["CSSPLIT"].split(",") for cs in midsv_control]

[cs.count("N") for cs in cssplits_sample]
[cs[1500] for cs in cssplits_sample]  # large
[cs[2200] for cs in cssplits_sample]  # all
[cs[2500] for cs in cssplits_sample]  # medium-large
[cs[3000] for cs in cssplits_sample]  # large

[cs[1901] for cs in cssplits_sample]  # large
[cs[1901] for cs in cssplits_control]  # large

cssplits_control = replace_both_ends_n(cssplits_control)
cssplits_sample = replace_both_ends_n(cssplits_sample)
mutation_score = make_score(cssplits_control, cssplits_sample)
scores_control = annotate_score(cssplits_control, mutation_score)
scores_sample = annotate_score(cssplits_sample, mutation_score)
labels = return_labels(scores_sample, scores_control)
# 1 Counter({3: 638, 2: 482, 1: 281, 16: 265, 0: 211, 18: 173, 13: 4, 8: 4, 6: 3, 11: 3, 15: 3, 5: 3, 4: 3, 7: 3, 9: 2, 12: 2, 17: 1, 10: 1, 14: 1})
Counter(labels)
# Counter({2: 500, 1: 476, 4: 460, 3: 9, 5: 5})
[cs[1500] for cs in scores_sample]  # large
[cs[2200] for cs in scores_sample]  # all
[cs[2500] for cs in scores_sample]  # medium-large
[cs[3000] for cs in scores_sample]  # large


def test_labels():
    for _ in range(1, 5):
        labels = GaussianMixture(n_components=10, random_state=0).fit_predict(X)
        print(Counter(labels))
    # Counter({8: 213, 4: 199, 0: 184, 5: 161, 7: 94, 9: 26, 6: 10, 2: 1, 3: 1, 1: 1})


test_labels()
test_labels()

for samp in scores_sample:
    sum(samp)


Counter(labels)


scaler = StandardScaler()
scores_scaler = scaler.fit_transform(scores_sample)
pca = PCA(n_components=3).fit(scores_scaler)
variance = pca.explained_variance_
X = pca.transform(scores_scaler) * variance

plt.scatter(X[:, 0], X[:, 1])
plt.xlabel("PC1")
plt.ylabel("PC2")
plt.show()

x = []
for mut in mutation_score:
    x.append(sum(mut.values()))

scores_sample[0]
mutation_score[1900]
x[1900]

labels = GaussianMixture(n_components=3, random_state=0).fit_predict(scores_sample)
Counter(labels)
labels = GaussianMixture(n_components=3, random_state=0).fit_predict(X)
Counter(labels)

plt.plot(scores_sample[0], label="medium")
plt.plot(scores_sample[1], label="small")
plt.plot(scores_sample[2], label="large")
plt.legend()
plt.show()
# =-----


"""
Controlが優先的に分かれていることが原因
-> Controlのなかで異常なサンプルを分ける
"""
i = 4
labels = GaussianMixture(n_components=i, random_state=0).fit_predict(X)
labels = labels.tolist()
labels_control = labels[len(scores_sample) :]
labels_sample = labels[: len(scores_sample)]
labels_merged = merge_clusters(labels_control, labels_sample)
labels_reorder = reorder_labels(labels_merged)
print(Counter(labels_reorder), Counter(labels_control))


[i for i, cs in enumerate(labels_control) if cs == 2]

X_control = reduce_dimention([], scores_control)
labels = GaussianMixture(n_components=20, random_state=0).fit_predict(X_control)
label_most = Counter(labels).most_common()[0][0]
scores_control_subset = [scores for label, scores in zip(labels, scores_control) if label == label_most][:1000]
len(scores_control_subset)


# 19 Counter({0: 623, 1: 454, 2: 428, 17: 290, 13: 205, 3: 186, 10: 68, 18: 58, 7: 4, 12: 4, 8: 4, 9: 3, 11: 3, 6: 3, 15: 2, 5: 2, 4: 2, 14: 1, 16: 1}) Counter({6: 454, 3: 428, 1: 290, 5: 186, 2: 58, 4: 34}) Counter({0: 620, 13: 203, 10: 68})
