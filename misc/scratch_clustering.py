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


def replace_repdel_to_n(count: list[dict], repeat_dels: set):
    count_replaced = deepcopy(count)
    for i, cnt in enumerate(count):
        if i not in repeat_dels:
            continue
        val_del = 0
        for cs, val in cnt.items():
            if cs.startswith("-"):
                val_del += val
                del count_replaced[i][cs]
        if "N" in count_replaced[i].keys():
            count_replaced[i]["N"] += val_del
        else:
            count_replaced[i].update({"N": val_del})
        return count_replaced


def compensate_n(count: list[dict]):
    count_replaced = deepcopy(count)
    for i, cnt in enumerate(count):
        if "N" not in cnt.keys():
            continue
        n_num = cnt["N"]
        elements = []
        probs = []
        coverage = sum(cnt.values()) - n_num
        for key, val in cnt.items():
            if key == "N":
                continue
            elements.append(key)
            probs.append(val / coverage)
        np.random.seed(1)
        samples = np.random.choice(a=elements, size=n_num, p=probs)
        for key in samples:
            count_replaced[i][key] += 1
        del count_replaced[i]["N"]
    return count_replaced


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

allele = "deletion"
sv = True
sequence = DICT_ALLELE[allele]
knockin_loci = KNOCKIN_LOCI[allele]

# Control
midsv_control = midsv.read_jsonl((Path(TEMPDIR, "midsv", f"{CONTROL_NAME}_{allele}.jsonl")))
cssplits_control = [cs["CSSPLIT"].split(",") for cs in midsv_control]
transpose_control = list(zip(*cssplits_control))
count_control = call_count(transpose_control)

# Sample
cssplits_sample = [cs["CSSPLIT"].split(",") for cs in classif_sample if cs["ALLELE"] == allele and cs["SV"] == sv]
transpose_sample = list(zip(*cssplits_sample))\
count_sample = call_count(transpose_sample)

# Find repetitive dels
repeat_dels = find_repetitive_dels(count_control, count_sample, sequence)
count_repdel_to_n_control = replace_repdel_to_n(count_control, repeat_dels)
count_repdel_to_n_sample = replace_repdel_to_n(count_sample, repeat_dels)

# Distribute N
count_compensated_control = compensate_n(count_repdel_to_n_control)
count_compensated_sample = compensate_n(count_repdel_to_n_sample)

# Mask common mutations
percent_control = call_percent(count_compensated_control)
percent_sample = call_percent(count_compensated_sample)
percent_subtraction = subtract_percentage(percent_control, percent_sample)
percent_discarded = discard_common_error(percent_subtraction, 0.5)
percent_discarded = discard_match(percent_discarded)

scores_control = annotate_scores(count_compensated_control, percent_discarded)
scores_sample = annotate_scores(count_compensated_sample, percent_discarded)

scores_control[0]
percent_control[1000]
percent_sample[1000]
percent_subtraction[1000]
percent_discarded[1000]

# qscores = [cs["QSCORE"].split(",") for cs in midsv_control]

# cssplits_replaced = replace_n_to_match(cssplits_control, sequence)
# cssplits_replaced = replace_lowquality_to_match(cssplits_replaced, qscores, sequence)
# transpose_control = transpose(cssplits_replaced)
# count_control = call_count(transpose_control)

# knockin_loci = find_knockin_loci(count_control, sequence)
# count_control = recount_knockin_loci(count_control, knockin_loci, sequence)

# # Sample
# cssplits_sample = [cs["CSSPLIT"].split(",") for cs in classif_sample if cs["ALLELE"] == allele and cs["SV"] == sv]
# transpose_sample = transpose(cssplits_sample)
# qscores = [cs["QSCORE"].split(",") for cs in classif_sample if cs["ALLELE"] == allele and cs["SV"] == sv]
# cssplits_replaced = replace_n_to_match(cssplits_sample, sequence)
# cssplits_replaced = replace_lowquality_to_match(cssplits_replaced, qscores, sequence)
# transpose_sample = transpose(cssplits_replaced)
# count_sample = call_count(transpose_sample)
