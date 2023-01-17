from __future__ import annotations

"""
課題
- Cables2 floxアレルにおいて、knock-in配列の補正がされていないため、knock-in配列上のシークエンスエラーによってアレルがクラスタリングされてしまっている
    - 1776番目の*GA変異
    - 2444番目の-T欠失
    - 2438の-C, 2456の-C, 2458の*TC
対策
- kcnok-in配列上のkmerをとり、コントロールのkmerと類似の配列を探す
    - FASTA_ALLELEのマッピングによりknock-in配列箇所を探す
    - knock-in配列を完全に含むように、周辺配列から1-2mer選択する
- 類似配列において、knock-in配列と変異プロファイルを比較し、変異の類似性が高い場合にはその変異を置き換える

メモ
- knock-in配列のスタート
    left-loxp: 1733
    right-loxp: 2428
結果
"""
from difflib import get_close_matches
import re
from collections import Counter, defaultdict
from copy import deepcopy
from pathlib import Path
import midsv
from scipy import stats
from scipy.spatial.distance import cosine
import midsv
from collections import defaultdict, Counter
from pathlib import Path
from src.DAJIN2.core import preprocess, classification, clustering, consensus, report
from itertools import groupby
from itertools import permutations
from collections import defaultdict
from pathlib import Path

from src.DAJIN2.core.preprocess import mappy_align

for fastas in permutations(Path("tests/data/preprocess/correct_knockin/").iterdir(), 2):
    print(fastas)
    ref, query = fastas
    ref_allele = ref.stem
    alignments = mappy_align.to_sam(ref, query, preset="splice")
    alignments = list(alignments)
    alignments = [a.split("\t") for a in alignments]
    alignments_midsv = midsv.transform(alignments, midsv=False, cssplit=True, qscore=False)[0]
    cssplits = alignments_midsv["CSSPLIT"].split(",")
    print(cssplits)
    mutations = dict()
    for i, cs in enumerate(cssplits):
        if cs.startswith("="):
            continue
        mutations.update({i: cs})

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
FASTA_ALLELES = preprocess.format_inputs.dictionize_allele(ALLELE)

allele = "flox"
sv = False
midsv_sample = midsv.read_jsonl(Path(TEMPDIR, "midsv", f"barcode31_splice_{allele}.jsonl"))
midsv_control = midsv.read_jsonl(Path(TEMPDIR, "midsv", f"barcode42_splice_{allele}.jsonl"))
midsv_sample[1]
midsv_control[3]

cssplits_sample = [m["CSSPLIT"].split(",") for m in midsv_sample]
cssplits_control = [m["CSSPLIT"].split(",") for m in midsv_control]
# cssplits_sample = []
# for samp in classif_sample:
#     if samp["ALLELE"] == allele and not samp["SV"]:
#         cssplits_sample.append(samp)

idx = 2444  # -A
idx = 1776  # *GA
# left_start = 1733
count_sample = defaultdict(int)
for cs in midsv_sample:
    x = cs["CSSPLIT"].split(",")[idx]
    count_sample[x] += 1

count_control = defaultdict(int)
for cs in midsv_control:
    x = cs["CSSPLIT"].split(",")[idx]
    count_control[x] += 1

count_sample
count_control

# cssplits_sample
# sequence = ["TAGCTAAT"]
# cssplits_sample = ["=T", "=A", "-A", "+G|C|T|=A", "A", "+C|-T", "=T"]

# cssplits_decompose = []
# for i, cssplits in enumerate(cssplits_sample):
#     if cssplits.startswith("-") or cssplits == "N":
#         continue


# knockin_kmer = ["AGCTA", "TAACT"]

###############################################################################
# functions
###############################################################################


def extract_knockin_loci(TEMPDIR) -> defaultdict[set]:
    """
    Extract differencial loci between alleles
        - The purpose is to lower match_score between very similar alleles such as point mutation.
    """
    fasta_alleles = list(Path(TEMPDIR, "fasta").iterdir())
    knockin_alleles = defaultdict(set)
    for comb in list(permutations(fasta_alleles, 2)):
        ref, query = comb
        ref_allele = ref.stem
        alignments = mappy_align.to_sam(ref, query, preset="splice")
        alignments = [a.split("\t") for a in alignments]
        alignments_midsv = midsv.transform(alignments, midsv=False, cssplit=True, qscore=False)[0]
        cssplits = alignments_midsv["CSSPLIT"].split(",")
        knockin = set()
        for i, cs in enumerate(cssplits):
            if cs == "N":
                knockin.add(i)
        knockin_alleles[ref_allele] = knockin
    return knockin_alleles


def split_sequence_in_5mer(sequence: str):
    sequence_kmer = dict()
    for i in range(len(sequence) - 5):
        sequence_kmer.update({sequence[i : i + 5]: i + 2})
    return sequence_kmer


def get_5mer_of_knockin_loci(sequence: str, knockin_loci: list):
    sequence_kmer = dict()
    for i in knockin_loci:
        sequence_kmer.update({sequence[i : i + 5]: i + 2})
    return sequence_kmer


def get_idx_of_similar_5mers(knockin_kmer: dict, sequence_kmer: dict) -> defaultdict(list):
    idx_of_similar_5mers = defaultdict(list)
    for kmer, locus in knockin_kmer.items():
        idxes = []
        for key in get_close_matches(kmer, sequence_kmer.keys(), n=100, cutoff=0.0):
            idx = sequence_kmer[key]
            if set(range(idx - 2, idx + 3)) & knockin_loci:
                continue
            idxes.append(idx)
        idx_of_similar_5mers[locus] = idxes
    return idx_of_similar_5mers


def count_indel_5mer(cssplits_transposed, indexes) -> defaultdict(dict):
    count_5mer = defaultdict(dict)
    for i in indexes:
        cssplits_5mer = cssplits_transposed[i - 2 : i + 3]
        count = {"ins": [1] * 5, "del": [1] * 5, "sub": [1] * 5}
        for j, cs in enumerate(cssplits_5mer):
            counter = Counter(cs)
            for key, cnt in counter.items():
                if key.startswith("=") or key == "N" or re.search(r"a|c|g|t|n", key):
                    continue
                if key.startswith("+"):
                    count["ins"][j] += cnt
                elif key.startswith("-"):
                    count["del"][j] += cnt
                elif key.startswith("*"):
                    count["sub"][j] += cnt
        count_5mer[i] = count
    return count_5mer


def replace_errors_to_match(cssplits_sample, sequence_errors, sequence):
    cssplits_replaced = []
    for cssplits in cssplits_sample:
        cssplits_copy = deepcopy(cssplits)
        for idx, error in sequence_errors.items():
            cssplits_5mer = cssplits_copy[idx - 2 : idx + 3]
            for j, mer in enumerate(cssplits_5mer):
                match_seq = "=" + sequence[idx - 2 + j]
                if "ins" in error and mer.startswith("+"):
                    cssplits_5mer[j] = match_seq
                if "del" in error and mer.startswith("-"):
                    cssplits_5mer[j] = match_seq
                if "del" in error and mer.startswith("*"):
                    cssplits_5mer[j] = match_seq
            cssplits_copy[idx - 2 : idx + 3] = cssplits_5mer
        cssplits_replaced.append(cssplits_copy)
    return cssplits_replaced


##########################################################
# main
# Todo 入力がFASTA_ALLELE、出力はmidsv.write_jsonl
##########################################################


def execute(TEMPDIR, FASTA_ALLELES):
    pass


knockin_alleles = extract_knockin_loci(TEMPDIR)
knockin_loci = knockin_alleles["flox"]
sequence = FASTA_ALLELES["flox"]
i = 1776

"""
コントロールとの類似配列において、エラーの相同性が高いものをシークエンスエラーとする
"""


knockin_kmer = get_5mer_of_knockin_loci(sequence, knockin_loci)
sequence_kmer = split_sequence_in_5mer(sequence)

idx_of_similar_5mers = get_idx_of_similar_5mers(knockin_kmer, sequence_kmer)

# i = 1776
# sequence[i - 2 : i + 3], i
# j = 1299
# sequence[j - 2 : j + 3]
# idx_of_similar_5mers[i]


count_5mer_similar_sequences = defaultdict(dict)
cssplits_transposed = [list(t) for t in zip(*cssplits_control)]
for i, indexes in idx_of_similar_5mers.items():
    count_5mer_similar_sequences[i] = count_indel_5mer(cssplits_transposed, indexes)

"""
Knock-in配列のエラープロファイル
"""
# i = 1776
transposed_cssplits = [list(t) for t in zip(*cssplits_sample)]
count_5mer_knockin = count_indel_5mer(transposed_cssplits, knockin_loci)


coverage_sample = len(midsv_sample)
coverage_control = len(midsv_control)

sequence_errors = defaultdict(set)

for i, count_knockin in count_5mer_knockin.items():
    samp_mutation = defaultdict(list)
    for mutation in ["ins", "del", "sub"]:
        samp_mutation[mutation] = [c / coverage_sample for c in count_knockin[mutation]]
    count_control = count_5mer_similar_sequences[i]
    for j, count in count_control.items():
        for mutation in ["ins", "del", "sub"]:
            samp = samp_mutation[mutation]
            cont = [c / coverage_control for c in count[mutation]]
            distance = 1 - cosine(samp, cont)
            _, pvalue = stats.ttest_ind(samp, cont, equal_var=False)
            if distance > 0.7 and pvalue > 0.01:
                sequence_errors[i].add(mutation)


###############################################################
# Knock-in配列中のシークエンスエラーを補正する
###############################################################

cssplits_replaced = replace_errors_to_match(cssplits_sample, sequence_errors, sequence)


# 確認

idx = 2444  # -A
idx = 1776  # *GA
# left_start = 1733
count_sample = defaultdict(int)
for cs in midsv_sample:
    x = cs["CSSPLIT"].split(",")[idx]
    count_sample[x] += 1

count_control = defaultdict(int)
for cs in midsv_control:
    x = cs["CSSPLIT"].split(",")[idx]
    count_control[x] += 1


count_replaced = defaultdict(int)
for cs in cssplits_replaced:
    x = cs[idx]
    count_replaced[x] += 1


count_sample
count_control
count_replaced


# i = 1776
# count_knockin = count_5mer_knockin[i]
# sequence_errors[i]
# count_control = count_control_5mer[73]
# count_control
# sorted(results, key=lambda x: -x[1])[:10]
# sorted(results, key=lambda x: -x[2])[:10]

###########################################################
# Knockin配列内において、エラーの相同性が高いものをシークエンスエラーとする
###########################################################

# cssplits = cssplits_sample
# transposed = [list(t) for t in zip(*cssplits)]
# count_knockin_5mer = defaultdict(dict)
# for i in knockin_loci:
#     cssplits_5mer = transposed[i - 2 : i + 3]
#     count = {"ins": [1] * 5, "del": [1] * 5, "sub": [1] * 5}
#     for j, cs in enumerate(cssplits_5mer):
#         counter = Counter(cs)
#         for key, cnt in counter.items():
#             if key.startswith("=") or key == "N" or re.search(r"a|c|g|t|n", key):
#                 continue
#             if key.startswith("+"):
#                 count["ins"][j] += cnt
#             elif key.startswith("-"):
#                 count["del"][j] += cnt
#             elif key.startswith("*"):
#                 count["sub"][j] += cnt
#     count_knockin_5mer[i] = count

# i = 2444
# results = []
# for i in knockin_loci:
#     for j in knockin_loci:
#         if j in range(i - 4, i + 5):
#             continue
#         x = count_5mer_knockin[i]
#         y = count_5mer_knockin[j]
#         op = "sub"
#         x = [c / coverage_sample for c in x[op]]
#         y = [c / coverage_sample for c in y[op]]
#         distance = 1 - cosine(x, y)
#         _, pvalue = stats.ttest_ind(x, y, equal_var=False)
#         # results.append((j, distance, pvalue))
#         if distance > 0.7 and pvalue > 0.01:
#             results.append((j, distance, pvalue))

# sorted(results, key=lambda x: -x[1])[:10]
# sorted(results, key=lambda x: -x[2])[:10]
# results
# sequence[i - 2 : i + 3]
# sequence[j - 2 : j + 3]
# sequence[i - 2 : j + 3]
# count_knockin_5mer[2492]
