from __future__ import annotations
from collections import Counter
import midsv
import numpy as np
import re
from pathlib import Path
from itertools import chain
from copy import deepcopy
import statsmodels.api as sm
import bisect
from collections import defaultdict
from pprint import pprint
from sklearn.cluster import MeanShift
from scipy.spatial.distance import cosine


# *Stx2 deletion: 2012-2739 (727 bases)
# tests/data/knockout/test_barcode25.fq.gz
# tests/data/knockout/test_barcode30.fq.gz
# tests/data/knockout/design_stx2.fa


def call_count(transpose_cssplits: list[list[str]]) -> list[dict[str:int]]:
    cssplit_counts = []
    for cssplit in transpose_cssplits:
        count = Counter(cssplit)
        count = dict(count)
        cssplit_counts.append(count)
    return cssplit_counts


def count_deletons_in_repaet(count_control, repeat_span):
    list_deletions = []
    count_deletions = []
    repeat_idx = 0
    repeat_start = repeat_span[repeat_idx][0]
    repeat_end = repeat_span[repeat_idx][1]
    for i, counts in enumerate(count_control):
        if repeat_start <= i < repeat_end:
            count_del = sum(val for key, val in counts.items() if key.startswith("-"))
            count_deletions.append(count_del)
        elif i == repeat_end:
            list_deletions.append(count_deletions)
            repeat_idx += 1
            if repeat_idx == len(repeat_span):
                break
            repeat_start = repeat_span[repeat_idx][0]
            repeat_end = repeat_span[repeat_idx][1]
            count_deletions = []
            if repeat_start == i:
                count_del = sum(val for key, val in counts.items() if key.startswith("-"))
                count_deletions.append(count_del)
    return list_deletions


def find_repetitive_dels(count_control, count_sample, sequence) -> set[int]:
    repeat_regrex = r"A{4,}|C{4,}|G{4,}|T{4,}|N{4,}"
    repeat_span = list(x.span() for x in re.finditer(repeat_regrex, sequence))
    control_repdels = count_deletons_in_repaet(count_control, repeat_span)
    sample_repdels = count_deletons_in_repaet(count_sample, repeat_span)
    cossims = [1 - cosine(c, s) for c, s in zip(control_repdels, sample_repdels)]
    high_cossims = [list(range(span[0], span[1])) for cossim, span in zip(cossims, repeat_span) if cossim > 0.95]
    return set(chain.from_iterable(high_cossims))


# def replace_repdel_to_n(count: list[dict], repeat_dels: set):
#     count_replaced = deepcopy(count)
#     for i, cnt in enumerate(count):
#         if i not in repeat_dels:
#             continue
#         val_del = 0
#         for cs, val in cnt.items():
#             if cs.startswith("-"):
#                 val_del += val
#                 del count_replaced[i][cs]
#         if "N" in count_replaced[i].keys():
#             count_replaced[i]["N"] += val_del
#         else:
#             count_replaced[i].update({"N": val_del})
#         return count_replaced


def replace_repdels_to_n(transposed_cssplits: list[list[str]], repeat_dels: set):
    replased_cssplits = deepcopy(transposed_cssplits)
    for i, cssplit in enumerate(replased_cssplits):
        if i not in repeat_dels:
            continue
        for j, cs in enumerate(cssplit):
            if cs.startswith("-"):
                replased_cssplits[i][j] = "N"
    return replased_cssplits


def sampling(cnt: Counter):
    elements = []
    probs = []
    coverage = sum(cnt.values()) - cnt["N"]
    for key, val in cnt.items():
        if key == "N":
            continue
        elements.append(key)
        probs.append(val / coverage)
    np.random.seed(1)
    samples = np.random.choice(a=elements, size=cnt["N"], p=probs)
    return samples


def replace_both_ends_n(transposed_cssplits: list[list[str]]):
    d_samples = defaultdict(iter)
    for i, cssplits in enumerate(transposed_cssplits):
        cnt = Counter(cssplits)
        if cnt["N"] == 0:  # No N
            continue
        if len(cnt) == 1 and cnt["N"] > 0:  # All N
            continue
        samples = sampling(cnt)
        d_samples[i] = iter(samples)
    cssplits_replaced = [list(cs) for cs in zip(*transposed_cssplits)]
    for i, cssplits in enumerate(cssplits_replaced):
        for j, cs in enumerate(cssplits):
            if cs != "N":
                break
            cssplits_replaced[i][j] = next(d_samples[j])
    for i, cssplits in enumerate(cssplits_replaced):
        cssplits = cssplits[::-1]
        for j, cs in enumerate(cssplits):
            if cs != "N":
                break
            cssplits_replaced[i][len(cssplits) - 1 - j] = next(d_samples[len(cssplits) - 1 - j])
    cssplits_replaced = cssplits_replaced[::-1]
    cssplits_replaced = [list(cs) for cs in zip(*cssplits_replaced)]
    return cssplits_replaced


# def compensate_n(count: list[dict]):
#     count_replaced = deepcopy(count)
#     for i, cnt in enumerate(count):
#         if "N" not in cnt.keys():
#             continue
#         n_num = cnt["N"]
#         elements = []
#         probs = []
#         coverage = sum(cnt.values()) - n_num
#         for key, val in cnt.items():
#             if key == "N":
#                 continue
#             elements.append(key)
#             probs.append(val / coverage)
#         np.random.seed(1)
#         samples = np.random.choice(a=elements, size=n_num, p=probs)
#         for key in samples:
#             count_replaced[i][key] += 1
#         del count_replaced[i]["N"]
#     return count_replaced


###############################################################################
# Discard common errors
###############################################################################


def call_percent(cssplit_counts: list[dict[str:int]]) -> list[dict[str:int]]:
    cssplit_percent = []
    coverage = sum(cssplit_counts[0].values())
    for counts in cssplit_counts:
        percent = {k: v / coverage * 100 for k, v in counts.items()}
        cssplit_percent.append(percent)
    return cssplit_percent


def subtract_percentage(percent_control, percent_sample) -> list[dict]:
    sample_subtracted = []
    for cont, samp in zip(percent_control, percent_sample):
        samp = Counter(samp)
        samp.subtract(Counter(cont))
        sample_subtracted.append(dict(samp))
    return sample_subtracted


def discard_common_error(sample_subtracted, threshold=0.5):
    sample_discarded = []
    for samp in sample_subtracted:
        remained = {k: v for k, v in samp.items() if v > threshold}
        sample_discarded.append(remained)
    return sample_discarded


def discard_match(sample_subtracted):
    sample_discarded = []
    for samp in sample_subtracted:
        remained = {k: v for k, v in samp.items() if not k.startswith("=")}
        sample_discarded.append(remained)
    return sample_discarded


###############################################################################
# annotate scores
###############################################################################


def annotate_scores(count_compensated, percent_discarded):
    scores = []
    for samp, mutation in zip(count_compensated, percent_discarded):
        score = []
        if not mutation:
            # # TODO 全インデックスにスコアを当てるとデバッグがし易い。本番では計算量低減のために以下の2行は削除する。
            # score = [0 for _ in range(len(samp))]
            # scores.append(score)
            continue
        for key, val in mutation.items():
            for s in samp:
                if s == key:
                    score.append(val)
                else:
                    score.append(0)
        scores.append(score)
    return list(zip(*scores))


###############################################################################

allele = "control"
sv = True
sequence = DICT_ALLELE[allele]
knockin_loci = KNOCKIN_LOCI[allele]

# Control
midsv_control = midsv.read_jsonl((Path(TEMPDIR, "midsv", f"{CONTROL_NAME}_{allele}.jsonl")))
cssplits_control = [cs["CSSPLIT"].split(",") for cs in midsv_control]
transpose_control = [list(cs) for cs in zip(*cssplits_control)]
count_control = call_count(transpose_control)

# Sample
cssplits_sample = [cs["CSSPLIT"].split(",") for cs in classif_sample if cs["ALLELE"] == allele and cs["SV"] == sv]
transpose_sample = [list(cs) for cs in zip(*cssplits_sample)]
count_sample = call_count(transpose_sample)

# Find repetitive dels
repeat_dels = find_repetitive_dels(count_control, count_sample, sequence)
transpose_repdels_control = replace_repdels_to_n(transpose_control, repeat_dels)
transpose_repdels_sample = replace_repdels_to_n(transpose_sample, repeat_dels)

# Distribute N
transpose_compensated_control = replace_both_ends_n(transpose_repdels_control)
transpose_compensated_sample = replace_both_ends_n(transpose_repdels_sample)

# Mask common mutations
count_compensated_control = call_count(transpose_compensated_control)
count_compensated_sample = call_count(transpose_compensated_sample)
percent_control = call_percent(count_compensated_control)
percent_sample = call_percent(count_compensated_sample)
percent_subtraction = subtract_percentage(percent_control, percent_sample)
percent_discarded = discard_common_error(percent_subtraction, 0.5)
percent_discarded = discard_match(percent_discarded)

scores_control = annotate_scores(transpose_compensated_control, percent_discarded)
scores_sample = annotate_scores(transpose_compensated_sample, percent_discarded)

len(scores_control)
len(scores_control[0])

# Clustering
scores_control_subset = scores_control[:1000]
scores = scores_sample + scores_control_subset
from time import time

time_start = time()
labels = MeanShift(n_jobs=THREADS).fit(scores).labels_
time_exec = time() - time_start
print(time_exec)
Counter(labels)

Counter(labels[: len(scores_sample)])

d = []
for s in scores_sample[:10]:
    count = sum(1 for ss in s if ss > 30)
    d.append(count)

d
