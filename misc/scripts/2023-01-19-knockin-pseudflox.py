"""
問題：最近の結果ではPseudo-floxアレルがなくなってしまったように見られる
"""

import re
from collections import Counter, defaultdict
from copy import deepcopy
from difflib import get_close_matches
from itertools import groupby, permutations
from pathlib import Path

import midsv
from scipy import stats
from scipy.spatial.distance import cosine

from src.DAJIN2.core import classification, clustering, consensus, preprocess, report
from src.DAJIN2.core.clustering import clustering
from src.DAJIN2.core.preprocess import mappy_align

# * flox insertion
SAMPLE, CONTROL, ALLELE, NAME, GENOME, DEBUG, THREADS = (
    "examples/flox-cables2/AyabeTask1/barcode31.fq.gz",
    "examples/flox-cables2/AyabeTask1/barcode42.fq.gz",
    "examples/flox-cables2/AyabeTask1/design_cables2.fa",
    "test-ayabe-task1",
    "mm10",
    True,
    14,
)

TEMPDIR = Path("DAJINResults", ".tempdir", NAME)
SAMPLE_NAME = preprocess.format_inputs.extract_basename(SAMPLE)
CONTROL_NAME = preprocess.format_inputs.extract_basename(CONTROL)
FASTA_ALLELES = preprocess.format_inputs.dictionize_allele(ALLELE)

if GENOME:
    UCSC_URL, GOLDENPATH_URL = preprocess.check_inputs.check_and_fetch_genome(GENOME)
    GENOME_COODINATES = preprocess.format_inputs.fetch_coodinate(GENOME, UCSC_URL, FASTA_ALLELES["control"])
    CHROME_SIZE = preprocess.format_inputs.fetch_chrom_size(GENOME_COODINATES["chr"], GENOME, GOLDENPATH_URL)

###############################################################################
allele = "flox"
sv = False
midsv_sample = midsv.read_jsonl(Path(TEMPDIR, "midsv", f"barcode31_splice_{allele}.jsonl"))
midsv_control = midsv.read_jsonl(Path(TEMPDIR, "midsv", f"barcode42_splice_{allele}.jsonl"))

both_flox = "f855fb9b57ee"
no_flox = "54661f075dc2"
left_flox = "7089ae32a686"
right_flox = "af9e73c7f917"

for m in midsv_sample:
    if left_flox in m["QNAME"]:
        print(m["QNAME"], m["CSSPLIT"])

"""bash
cat DAJINResults/.tempdir/test-ayabe-task1/sam/barcode31_splice_flox.sam |
    grep 7089ae32a686 | grep -c "~"
cat DAJINResults/.tempdir/test-ayabe-task1/sam/barcode31_splice_flox.sam |
    grep af9e73c7f917 | grep -c "~"
"""

classif_sample = classification.classify_alleles(TEMPDIR, SAMPLE_NAME)

for classif in classif_sample:
    classif["SV"] = classification.detect_sv(classif["CSSPLIT"], threshold=50)

for c in classif_sample:
    if left_flox in c["QNAME"]:
        print(c)  # left_floxはflox, sv: false

for c in classif_sample:
    if right_flox in c["QNAME"]:
        print(c)  # right_floxはcontrol, sv:false

clust_sample = clustering.add_labels(classif_sample, TEMPDIR, CONTROL_NAME, FASTA_ALLELES, THREADS)
clust_sample = clustering.add_readnum(clust_sample)
clust_sample = clustering.add_percent(clust_sample)
clust_sample = clustering.update_labels(clust_sample)

for c in clust_sample:
    if both_flox in c["QNAME"]:
        print(c)  # both_floxはLabel1


for c in clust_sample:
    if left_flox in c["QNAME"]:
        print(c)  # left_floxはLabel1

for c in clust_sample:
    if right_flox in c["QNAME"]:
        print(c)  # right_floxはlabel15

for c in clust_sample:
    if no_flox in c["QNAME"]:
        print(c)  # no_floxはlabel4

###########################################################
# Classification
###########################################################
from src.DAJIN2.core.clustering.replace_both_ends_n import replace_both_ends_n
from src.DAJIN2.core.clustering.make_score import make_score
from src.DAJIN2.core.clustering.annotate_score import annotate_score
from src.DAJIN2.core.clustering.reorder_labels import reorder_labels
from src.DAJIN2.core.clustering.return_labels import return_labels

paths_midsv = list(Path(TEMPDIR, "midsv").glob(f"{CONTROL_NAME}_splice_*"))
cssplits_control_by_alleles = defaultdict(list)
for path_midsv in paths_midsv:
    midsv_control = midsv.read_jsonl(path_midsv)
    cssplits = [cs["CSSPLIT"].split(",") for cs in midsv_control]
    allele = path_midsv.stem.split("_")[-1]
    cssplits_control_by_alleles[allele] = cssplits


allele = "flox"
sv = False

cssplits_control = cssplits_control_by_alleles[allele]
cssplits_sample = [cs["CSSPLIT"].split(",") for cs in classif_sample if cs["ALLELE"] == allele and cs["SV"] == sv]
cssplits_control = replace_both_ends_n(cssplits_control)
cssplits_sample = replace_both_ends_n(cssplits_sample)
mutation_score = make_score(cssplits_control, cssplits_sample)
scores_control = annotate_score(cssplits_control, mutation_score)
scores_sample = annotate_score(cssplits_sample, mutation_score)
labels = return_labels(scores_sample, scores_control)
Counter(labels)

sequence = FASTA_ALLELES[allele]
len(sequence)
len(mutation_score)

samp = []
for i, cs in enumerate(classif_sample):
    if cs["ALLELE"] == allele and cs["SV"] == sv:
        samp.append(cs)

for i, cs in enumerate(samp):
    if left_flox in cs["QNAME"]:
        print(i, "left")  # 212
    if both_flox in cs["QNAME"]:
        print(i, "both")  # 464

scores_sample[212]
scores_sample[464]
labels[212]
labels[464]
mutation_score[1740]

for i, l in enumerate(labels_sample):
    if l == 2:
        print(scores_sample[i])
        break

samp[i]

from matplotlib import pyplot as plt

sum(scores_control_subset[0])
for l, s in zip(labels_control, scores_control_subset):
    if l == 3:
        print(s)
        sum(s)
        break

fig, axs = plt.subplots(4)
axs[0].plot(scores_control_subset[0])
axs[1].plot(s)
axs[2].plot(scores_sample[212])
axs[3].plot(scores_sample[464])
plt.show(block=False)

s.index(-12.717271641883617)
cnt = Counter()
for cs in cssplits_control:
    cnt += Counter([cs[86]])

cnt
mutation_score[86]
scores_control_subset[86]
# for c in labels:
#     if both_flox in c["QNAME"]:
#         print(c)  # both_floxはLabel1

# for c in clust_sample:
#     if left_flox in c["QNAME"]:
#         print(c)  # left_floxはLabel1


###########################################################
# cssplits_controlになぜかpseudo-floxがある？
# -> pseudo-floxではなく、spliceによってspliceとなるか欠失となるかの差異のせいだった
###########################################################

cnt = defaultdict(int)
for cs in cssplits_control:
    num = len(list(re.finditer(r"N{10,}", "".join(cssplits_control[0]))))
    cnt[num] += 1

cnt

"""
cssplits_controlのなかにはpseudo-floxはない。
-> コードの中でcontrolとsampleが混ざってしまっている？
擬似的なscoers_sampleとscores_controlを作って検討
"""

# import random

# scores_sample = []
# for _ in range(500):
#     scores = []
#     for _ in range(5):
#         scores.append(random.randrange(1, 9))
#     scores_sample.append(scores)

# scores_control = []
# for _ in range(1000):
#     scores = []
#     for _ in range(5):
#         scores.append(random.randrange(1, 9))
#     scores_control.append(scores)

# labels = return_labels(scores_sample, scores_control)
# Counter(labels)

"""
scores_controlのなかにすでにpseud-floxが混じっている？
"""

cssplits_control = cssplits_control_by_alleles[allele]
cssplits_sample = [cs["CSSPLIT"].split(",") for cs in classif_sample if cs["ALLELE"] == allele and cs["SV"] == sv]
cssplits_control = replace_both_ends_n(cssplits_control)
cssplits_sample = replace_both_ends_n(cssplits_sample)
mutation_score = make_score(cssplits_control, cssplits_sample)
scores_control = annotate_score(cssplits_control, mutation_score)
scores_sample = annotate_score(cssplits_sample, mutation_score)

cnt = defaultdict(int)
for i, cs in enumerate(cssplits_control):
    num = len(list(re.finditer(r"N{10,}", "".join(cs))))
    cnt[num] += 1
    if num == 1:
        print(i)

cnt
cssplits_control[0]
cssplits_control[18]
midsv_control[18]


scores_mark = []
for score in scores_control:
    s = ["N" if s < -50 else "@" for s in score]
    scores_mark.append(s)

cnt = defaultdict(int)
for cs in scores_mark:
    num = len(list(re.finditer(r"N{10,}", "".join(cssplits_control[0]))))
    cnt[num] += 1

cnt


path_midsv = "DAJINResults/.tempdir/test-ayabe-task1/midsv/barcode42_splice_flox.jsonl"
midsv_control = midsv.read_jsonl(path_midsv)
cssplits = [cs["CSSPLIT"].split(",") for cs in midsv_control]
allele = path_midsv.stem.split("_")[-1]
cssplits_control_by_alleles[allele] = cssplits

midsv_control[18]

"""bash
cat examples/flox-cables2/AyabeTask1/design_cables2.fa | grep flox -A 1 > tmp_flox.fa
zcat examples/flox-cables2/AyabeTask1/barcode42.fq.gz | grep 055e927e-f4dc-44f9-aa17-862f0e3e68b0 -A 3 |
    minimap2 -ax map-ont --splice --cs=long tmp_flox.fa - > tmp_flox_splice.sam
    samtools sort tmp_flox_splice.sam > tmp_flox_splice.bam
    samtools index tmp_flox_splice.bam

cat DAJINResults/.tempdir/test-ayabe-task1/sam/barcode42_splice_flox.sam | grep 862f0e3e68b0 | tr "ACGT" "." | grep -e "~"
cat DAJINResults/.tempdir/test-ayabe-task1/midsv/barcode42_splice_flox.jsonl | grep 862f0e3e68b0 | grep -
"""


"""
僅かな欠失塩基数の差によってleft_loxpがスプライス（N）となるか欠失(-)となるか分かれており、これがcontrolのクラスタリング分けに影響されている
およそ60リード(968リード中) が欠失として判定されている

"""
scores_control[0][1750]
scores_control[18][1750]

for i, l in enumerate(labels):
    if l == 9:
        print(i)

scores_control[36]

fig, axs = plt.subplots(3)
axs[0].plot(scores_control[36])
axs[1].plot(scores_control[0])
axs[2].plot(scores_control[18])
plt.show(block=False)

labels[0]
labels[18]
labels[36]


# targets = [0, 18, 36]
# fig, axs = plt.subplots(len(targets))
# for i, t in enumerate(targets):
#     axs[i].plot(scores_control[t])

# plt.show(block=False)


from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA


scores = scores_control
n_components = min(20, len(scores[0]))
# scaler = StandardScaler()
# scores = scaler.fit_transform(scores)
pca = PCA(n_components=n_components).fit(scores)
# variance = pca.explained_variance_


x = pca.transform(scores)
x = pca.transform(scores)
targets = [0, 18, 36]
fig, axs = plt.subplots(len(targets))
for i, t in enumerate(targets):
    axs[i].plot(x[t])

plt.show(block=False)


labels = GaussianMixture(n_components=5, random_state=1).fit_predict(x)
Counter(labels)

labels[0]
labels[18]
labels[36]


targets = [0, 3]
fig, axs = plt.subplots(len(targets))
for i, t in enumerate(targets):
    axs[i].plot(scores_control_subset[t])

plt.show(block=False)

scores[len(scores_sample) :][0]

scores[212]
scores[464]
scores[len(scores_sample) + 1]
scores[len(scores_sample) + 4]

targets = [212, 464, len(scores_sample) + 1, len(scores_sample) + 4]
fig, axs = plt.subplots(len(targets))
for i, t in enumerate(targets):
    axs[i].plot(scores[t])

plt.show(block=False)

targets = [212, 464, len(scores_sample) + 1, len(scores_sample) + 4]
fig, axs = plt.subplots(len(targets))
for i, t in enumerate(targets):
    axs[i].plot(X[t])

plt.show(block=False)

labels[212]
labels[464]
mutation_score[1740]

"""
2023-01-20 次回への課題
- PCAの前処理によってかなり結果が変わる。
    - とくにスケーリングするか否かが重要であり、この挙動に対する理解が必要
    - いままではvarianceを掛けていたが、これではほぼPC1の情報しか残らないので掛けないようにした
        - 細かいノイズを拾ってしまう可能性があるので、ほかのアルビノサンプルなどで検討する必要がある
- controlにおけるGaussianMixtureのn_conponentsも重要。20だと分かれすぎるきらいがあり、いまは5としている
"""
