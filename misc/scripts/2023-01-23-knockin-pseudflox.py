"""
問題：ControlにアライメントされたPseudo-floxアレルがクラスタリングされていない
control="0c83268d26d0"
right_loxp1="6fe3a1bd5b88"
right_loxp2="22b9661586bf"
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
# midsv_sample = midsv.read_jsonl(Path(TEMPDIR, "midsv", f"barcode31_splice_{allele}.jsonl"))
# midsv_control = midsv.read_jsonl(Path(TEMPDIR, "midsv", f"barcode42_splice_{allele}.jsonl"))

classif_sample = classification.classify_alleles(TEMPDIR, SAMPLE_NAME)

for classif in classif_sample:
    classif["SV"] = classification.detect_sv(classif["CSSPLIT"], threshold=50)


allele = "control"
sv = False
intact = "0c83268d26d0"
right_loxp1 = "6fe3a1bd5b88"
right_loxp2 = "ee0f82c5d57"

[c for c in classif_sample if intact in c["QNAME"]]
[c for c in classif_sample if right_loxp1 in c["QNAME"]]
[c for c in classif_sample if right_loxp2 in c["QNAME"]]

"""bash
"""

clust_sample = clustering.add_labels(classif_sample, TEMPDIR, CONTROL_NAME, FASTA_ALLELES, THREADS)
clust_sample = clustering.add_readnum(clust_sample)
clust_sample = clustering.add_percent(clust_sample)
clust_sample = clustering.update_labels(clust_sample)

[c for c in clust_sample if intact in c["QNAME"]]  # Label2
[c for c in clust_sample if right_loxp1 in c["QNAME"]]  # Label2...
[c for c in clust_sample if right_loxp2 in c["QNAME"]]  # Label2...

###########################################################
# Classification
###########################################################
from DAJIN2.core.clustering.preprocess import replace_both_ends_n
from DAJIN2.core.clustering.preprocess import compress_insertion
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


allele = "control"
sv = False

cssplits_control = cssplits_control_by_alleles[allele]
cssplits_sample = [cs["CSSPLIT"].split(",") for cs in classif_sample if cs["ALLELE"] == allele and cs["SV"] == sv]
cssplits_control = replace_both_ends_n(cssplits_control)
cssplits_sample = replace_both_ends_n(cssplits_sample)
cssplits_control = compress_insertion(cssplits_control)
cssplits_sample = compress_insertion(cssplits_sample)
mutation_score = make_score(cssplits_control, cssplits_sample)
scores_control = annotate_score(cssplits_control, mutation_score)
scores_sample = annotate_score(cssplits_sample, mutation_score)
labels = return_labels(scores_sample, scores_control)
Counter(labels)


samp = []
for i, cs in enumerate(classif_sample):
    if cs["ALLELE"] == allele and cs["SV"] == sv:
        samp.append(cs)

for i, cs in enumerate(samp):
    if intact in cs["QNAME"]:
        print(i, "intact")  # 344
    if right_loxp1 in cs["QNAME"]:
        print(i, "right_loxp1")  # 10

samp[334]["CSSPLIT"].count("+")


count_loxp = 0
indexes_loxp = []
for i, s in enumerate(samp):
    if re.search(r"(\+[ACGT]\|){10}", s["CSSPLIT"]):
        count_loxp += 1
        indexes_loxp.append(i)

len(samp) - count_loxp, count_loxp

scores_sample[344]
scores_sample[10]
samp[11]["CSSPLIT"].split(",")[1730:1740]
cssplits_sample[11][1730:1740]
labels[344]
labels[10]

for i, s in enumerate(samp[10]["CSSPLIT"].split(",")):
    if "+C|+A|+A|+A|+C|+A|+T|" in s:
        print(i)  # 1732

percent_sample[1732]
percent_control[1732]
mutation_score[1732]
scores_sample[10][1731:1734]
scores_control[10][1731:1734]


"""
挿入塩基中のミスマッチによって細切れに分かれてしまっている
→ CSSPLIT中の挿入を、挿入塩基ではなく「挿入塩基数」に変更する
"""

cssplits_control = cssplits_control_by_alleles[allele]
cssplits_sample = [cs["CSSPLIT"].split(",") for cs in classif_sample if cs["ALLELE"] == allele and cs["SV"] == sv]
cssplits_control = replace_both_ends_n(cssplits_control)
cssplits_sample = replace_both_ends_n(cssplits_sample)
cssplits_control = compress_insertion(cssplits_control)
cssplits_sample = compress_insertion(cssplits_sample)
mutation_score = make_score(cssplits_control, cssplits_sample)
scores_control = annotate_score(cssplits_control, mutation_score)
scores_sample = annotate_score(cssplits_sample, mutation_score)
labels = return_labels(scores_sample, scores_control)
Counter(labels)


cssplits_sample[10][1732]

percent_sample[1732]
percent_control[1732]
percent_subtraction[1732]
percent_discarded[1732]
mutation_score[1732]
scores_sample[10][1732]
scores_sample[344][1732]
scores_sample[10][1732]
labels[344]
labels[10]


from matplotlib import pyplot as plt

fig, axs = plt.subplots(2)
axs[0].plot(scores_control[0])
axs[1].plot(scores_control[3])
plt.show(block=False)

scores_control[3][310:313]
sorted(mutation_score[313].items(), key=lambda x: -abs(x[1]))

sorted(percent_sample[313].items(), key=lambda x: -abs(x[1]))
sorted(percent_control[313].items(), key=lambda x: -abs(x[1]))

FASTA_ALLELES[allele][313:340]
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

"""
# 2023-01-26 進捗まとめ
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
