"""
# 2023-01-26 進捗まとめ 再掲
## 進捗箇条書き
- 2023-01-19に行った、`discard_common_error`のなかでsample - controlに絶対値を設けた
    - これによりcontrolにより多く存在する変異が原因でクラスタリングされることになった
    - そのため、sampleよりもむしろcontrolのほうが優先的にクラスタリングされるようになり、頑強性が小さくなってしまった…
    - # * Insertionを圧縮するアイディアはそのままでよい
        - # ? 一方でInsertionのなかにある変異を同定することができないので、これには別の手法が必要？
## 今後の課題
- # TODO floxアレルのpeudo-loxpからやり直す
- # TODO `discard_common_error`のなかにあるsample - controlの絶対値を外す
    - これによりcontrolに多くある変異を無視する
        - 我々が興味があるのはsample特異的な変異であり、controlの変異は興味がないため
    - しかし、これではfloxアレルのpseudo-loxpが同定できない
        - pseudo-loxpがないところはcontrolとおなじ`N`であるため、引き算したら消えてしまうため
        - # TODO loxpがあるところは、controlから引き算をしない（無視する）という形に置き換える
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
both_flox = "f855fb9b57ee"
left_flox = "7089ae32a686"

classif_sample = classification.classify_alleles(TEMPDIR, SAMPLE_NAME)

for classif in classif_sample:
    classif["SV"] = classification.detect_sv(classif["CSSPLIT"], threshold=50)

for c in classif_sample:
    if both_flox in c["QNAME"]:
        print(c)  # both_floxはallele: flox, sv: false

for c in classif_sample:
    if left_flox in c["QNAME"]:
        print(c)  # left_floxはallele: flox, sv: false


###########################################################
# Classification
###########################################################
from src.DAJIN2.core.preprocess.correct_knockin import extract_knockin_loci
from src.DAJIN2.core.clustering.preprocess import replace_both_ends_n
from src.DAJIN2.core.clustering.make_score import make_score
from src.DAJIN2.core.clustering.annotate_score import annotate_score
from src.DAJIN2.core.clustering.return_labels import return_labels

paths_midsv = list(Path(TEMPDIR, "midsv").glob(f"{CONTROL_NAME}_splice_*"))
cssplits_control_by_alleles = defaultdict(list)
for path_midsv in paths_midsv:
    midsv_control = midsv.read_jsonl(path_midsv)
    cssplits = [cs["CSSPLIT"].split(",") for cs in midsv_control]
    allele = path_midsv.stem.split("_")[-1]
    cssplits_control_by_alleles[allele] = cssplits

"""
Knockin_alleleを無視するために、knockin配列領域を入手する
preprocessにある`extract_knockin_loci`を再利用する
"""

knockin_alleles = extract_knockin_loci(TEMPDIR)

allele = "flox"
sv = False

knockin_loci = knockin_alleles[allele]
cssplits_control = cssplits_control_by_alleles[allele]
cssplits_sample = [cs["CSSPLIT"].split(",") for cs in classif_sample if cs["ALLELE"] == allele and cs["SV"] == sv]
cssplits_control = replace_both_ends_n(cssplits_control)
cssplits_sample = replace_both_ends_n(cssplits_sample)
"""
mutation_scoreの計算にknockin配列領域の情報を利用する
"""
mutation_score = make_score(cssplits_control, cssplits_sample, knockin_loci)
scores_control = annotate_score(cssplits_control, mutation_score)
scores_sample = annotate_score(cssplits_sample, mutation_score)
labels = return_labels(scores_sample, scores_control)
Counter(labels)

samp = []
for i, cs in enumerate(classif_sample):
    if cs["ALLELE"] == allele and cs["SV"] == sv:
        samp.append(cs)

for i, cs in enumerate(samp):
    if left_flox in cs["QNAME"]:
        print(i, "left")  # 212
    if both_flox in cs["QNAME"]:
        print(i, "both")  # 464

scores_sample[0]
scores_sample[212]
scores_sample[464]
labels[212]
labels[464]
mutation_score[1740]
scores_sample[212][1740]
scores_sample[464][1740]
cssplits_sample[212][1740]  # pseudo-leftloxp
cssplits_sample[464][1740]  # flox
percent_sample[1740]
percent_control[1740]


for i, m in enumerate(mutation_score):
    if any(True for v in m.values() if v > 5):
        print(i, m)


import urllib.request


def print_bed(sequence, idx, length=50):
    query_seq = sequence[idx : idx + length]
    query = f"https://gggenome.dbcls.jp/ja/mm10/{query_seq}.bed"
    with urllib.request.urlopen(query) as f:
        html = f.read().decode("utf-8")
    print(query_seq)
    print(*html.split("\n")[1].split("\t"))


print_bed(sequence, 2507)
x = [c[1] for c in point_coodinates]
y = [c[2] for c in point_coodinates]

from matplotlib import pyplot as plt

plt.plot(x, y, "o")
plt.show(block=False)

"""
Controlのクラスタ数が１を保ちつつ、サンプルが最大いくつに分かれるか、という設定のほうがノイズが取れなくて良いかも？
"""

Counter(labels_sample)
zero = labels_sample.index(0)
one = labels_sample.index(1)
two = labels_sample.index(2)
# [i for i,l in enumerate(labels_sample) if l == 2]
scores_sample[zero]
scores_sample[one]
scores_sample[two]


def return_consensus_label(cssplits, scores, labels, target_label):
    consensus_sequence = [dict()] * len(scores[0])
    consensus_scores = [0] * len(scores[0])
    cssplits_target = [cs for cs, label in zip(cssplits, labels) if label == target_label]
    scores_target = [s for s, label in zip(scores, labels) if label == target_label]
    count = len(cssplits_target)
    cssplits_t = list(zip(*cssplits_target))
    scores_t = list(zip(*scores_target))
    for i in range(len(cssplits_t)):
        consensus_sequence[i] = dict(Counter(cssplits_t[i]))
        consensus_scores[i] = sum(scores_t[i]) / count
    consensus_result = [(i, c, s) for i, (c, s) in enumerate(zip(consensus_sequence, consensus_scores)) if s > 1]
    return consensus_result


Counter(labels_sample)
scores_zero = return_consensus_label(scores_sample, labels_sample, 0)
scores_one = return_consensus_label(scores_sample, labels_sample, 1)
scores_two = return_consensus_label(scores_sample, labels_sample, 2)

scores_zero
scores_one
scores_two

###############################
# Controlの異常検知を試す
###############################

from sklearn.svm import OneClassSVM

X_control = reduce_dimention([], scores_control)
labels = OneClassSVM(gamma="auto").fit_predict(X_control)

cssplits, scores, labels, target_label = cssplits_control, scores_control, labels, 1
cont_normal = return_consensus_label(cssplits_control, scores_control, labels, 0)
cont_abnormal = return_consensus_label(cssplits_control, scores_control, labels, 1)

cont_normal
cont_abnormal
[i for i, l in enumerate(labels) if l == 1]
cssplits_control[19][1760:1780]
["" if cs.startswith("=") else cs for cs in cssplits_control[45]]
cssplits_control[0][1740]
cssplits_control[1][1740]
scores_control[0][1740]
scores_control[1][1740]
