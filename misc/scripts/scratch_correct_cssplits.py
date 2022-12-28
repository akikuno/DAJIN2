from __future__ import annotations
from collections import Counter
import re
from itertools import chain
from copy import deepcopy
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy import stats
from scipy.spatial.distance import cosine

import midsv

"""
- 5-merでコントロールと比較して、コントロールにもあるプロファイルの場合はマッチで補正する
    - 断片のリードで、端から続くNは無視
    - 補正する前にすでにすべてがマッチなら即continue
"""


def extract_indexes_with_both_ends_not_N(cssplits: list[list[str]]) -> list[tuple[int, int]]:
    indexes = []
    for cssplit in cssplits:
        cssplit = ",".join(cssplit)
        n_prefix = re.search(r"^(N,)+", cssplit)
        left_idx = n_prefix.end() if n_prefix else 0
        n_suffix = re.search(r"(,N)+$", cssplit)
        right_idx = n_suffix.start() if n_suffix else len(cssplit)
        # output index of splitted cssplit
        left_idx = cssplit[:left_idx].count(",")
        right_idx = cssplit.count(",") - cssplit[right_idx:].count(",")
        indexes.append((left_idx, right_idx))
    return indexes


def call_count(cssplits: list[list[str]], indexes: list[tuple(int, int)]) -> dict[dict[str, int]]:
    """Count cssplits within 3-mer range at mutation (or sequence error) loci.
    """
    count_kmer = defaultdict(Counter)
    for cssplit, idx in zip(cssplits, indexes):
        left_idx, right_idx = idx
        for i in range(left_idx + 1, right_idx):
            if cssplit[i].startswith("="):
                continue
            kmer = ",".join([cssplit[i - 1], cssplit[i], cssplit[i + 1]])
            count_kmer[i] += Counter([kmer])
    counts = {i: dict(count_kmer[i]) for i in count_kmer.keys()}
    return counts


def call_percentage(cssplits: list[list[str]], counts: dict[dict[str, int]]) -> dict[dict[str, float]]:
    coverage = len(cssplits)
    percents = deepcopy(counts)
    for i, c in counts.items():
        for kmer, count in c.items():
            percents[i][kmer] = count / coverage * 100
    return percents


def subtract_percentage(percent_sample: dict, percent_control: dict) -> dict[dict[str, float]]:
    percent_subtracted = deepcopy(percent_sample)
    for i, samp in percent_sample.items():
        cont = percent_control.get(i)
        for kmer_samp, per_samp in samp.items():
            if cont.get(kmer_samp):
                percent_subtracted[i][kmer_samp] = per_samp - cont[kmer_samp]
    return percent_subtracted


def select_candidate_mutation(percent_subtracted: dict, threshold: float = 0.5) -> dict[dict[str, float]]:
    candidate_mutation = dict()
    for i, samp in percent_subtracted.items():
        mutation = {k for k, v in samp.items() if v > threshold}
        candidate_mutation.update({i: mutation})
    return candidate_mutation


def update_cssplits(cssplits: list, sequence: str, candidate_mutation: dict) -> list[list[str]]:
    cssplits_update = deepcopy(cssplits)
    for j, mutation in candidate_mutation.items():
        for i, cssplit in enumerate(cssplits):
            kmer = ",".join([cssplit[j - 1], cssplit[j], cssplit[j + 1]])
            if kmer in mutation:
                continue
            else:
                cssplits_update[i][j] = "=" + sequence[j]
    return cssplits_update


###############################################################################
# Scratch... 5mers
###############################################################################

allele = "control"
sequence_length = len(FASTA_ALLELES[allele])
midsv_sample = midsv.read_jsonl("DAJINResults/.tempdir/tyr-pm/midsv/barcode31_splice_control.jsonl")
midsv_control = midsv.read_jsonl("DAJINResults/.tempdir/tyr-pm/midsv/barcode32_splice_control.jsonl")

cssplits_sample = [cs["CSSPLIT"].split(",") for cs in midsv_sample]
cssplits_control = [cs["CSSPLIT"].split(",") for cs in midsv_control]

sequence_length = 10
test_sample = []
for cs in cssplits_sample:
    test_sample.append(cs[51 - 9 : 52])

test_control = []
for cs in cssplits_control:
    test_control.append(cs[51 - 9 : 52])

"""
- Ins, Del, Sub, Invについてカウントする
"""
num_subset = sequence_length % 5
left_idx = 0
right_idx = sequence_length
if num_subset == 1:
    left_idx += 1
elif num_subset == 2:
    left_idx += 1
    right_idx -= 1
elif num_subset == 3:
    left_idx += 2
    right_idx -= 1
elif num_subset == 4:
    left_idx += 2
    right_idx -= 2

for t in test_sample:
    for i in range(left_idx, right_idx, 5):
        t[i : i + 5]

for i in range(left_idx, right_idx, 5):
    print(i, t[i])


def count_indels_5mer(cssplits: list[list[str]], left_idx: int, right_idx: int) -> list[dict]:
    transposed = [list(t) for t in zip(*cssplits)]
    count_indels_5mer = []
    for i in range(left_idx, right_idx, 5):
        count = {"ins": [1] * 5, "del": [1] * 5, "sub": [1] * 5}
        cssplits_5mer = transposed[i : i + 5]
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
        count_indels_5mer.append(count)
    return count_indels_5mer


count_5mer_sample = count_indels_5mer(test_sample)
count_5mer_control = count_indels_5mer(test_control)


def extractr_sequence_errors(count_5mer_sample, count_5mer_control):
    sequence_errors = [set() for _ in range(len(count_5mer_sample))]
    for i in range(len(sequence_errors)):
        for ids in ["ins", "del", "sub"]:
            x = count_5mer_sample[i][ids]
            y = count_5mer_control[i][ids]
            distance = 1 - cosine(x, y)
            _, pvalue = stats.ttest_ind(x, y, equal_var=False)
            if distance > 0.95 and pvalue > 0.05:
                sequence_errors[i].add(ids)
    return sequence_errors


sequence_errors = extractr_sequence_errors(count_5mer_sample, count_5mer_control)


def replace_errors_to_atmark(cssplits_sample, sequence_errors, left_idx, right_idx):
    cssplits_replaced = []
    for samp in cssplits_sample:
        samp_replaced = deepcopy(samp)
        for idx_error, idx_5mer in enumerate(range(left_idx, right_idx, 5)):
            samp_5mer = samp[idx_5mer : idx_5mer + 5]
            error = sequence_errors[idx_error]
            if "ins" in error:
                samp_5mer = ["@" if cs.startswith("+") else cs for cs in samp_5mer]
            if "del" in error:
                samp_5mer = ["@" if cs.startswith("-") else cs for cs in samp_5mer]
            if "sub" in error:
                samp_5mer = ["@" if cs.startswith("*") else cs for cs in samp_5mer]
            samp_replaced[idx_5mer : idx_5mer + 5] = samp_5mer
        cssplits_replaced.append(samp_replaced)
    return cssplits_replaced


cssplits_replaced = replace_errors_to_atmark(test_sample, sequence_errors, left_idx, right_idx)

for cs in cssplits_replaced:
    if any(c.startswith("-") for c in cs[0:5]):
        print(cs)


def correct_cssplits():
    """
    - samplingによって@をNとinv以外のMIDSに置き換える。
    - このsannplingも、5merにする？
    """
    pass


# d = defaultdict(int)
# for cs in cssplits_sample:
#     d[cs[51]] += 1

# d
# cssplits_sample[0][51]
