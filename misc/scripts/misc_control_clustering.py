# controlもクラスタリングして、controlにもあるクラスタは削除する

from __future__ import annotations
from collections import Counter
import midsv
import numpy as np
from sklearn.cluster import OPTICS
import re
from pathlib import Path
from itertools import chain
from copy import deepcopy
import statsmodels.api as sm
import bisect
from collections import defaultdict
from pprint import pprint

np.random.seed(1)

###############################################################################
# Replace N and Low Quality
###############################################################################


def replace_n_to_match(cssplits: list[list[str]], sequence: str) -> list[str]:
    for idx, cssplit in enumerate(cssplits):
        # Replace N to @ at the left ends
        for i, cs in enumerate(cssplit):
            if cs != "N":
                break
            cssplit[i] = "=" + sequence[i]
        # Replace N to @ at the right ends
        cssplit = cssplit[::-1]
        for i, cs in enumerate(cssplit):
            if cs != "N":
                break
            cssplit[i] = "=" + sequence[::-1][i]
        cssplits[idx] = cssplit[::-1]
    return cssplits


# cssplits = [["N", "N", "=A", "N", "N"]]
# sequence = "GCATG"
# replace_n_to_match(cssplits, sequence)


def replace_lowquality_to_match(cssplits: list[list[str]], qscores: list[list[int]], sequence: str) -> list[str]:
    for i, qscore in enumerate(qscores):
        for j, qs in enumerate(qscore):
            if "|" in qs:
                cs_ins = cssplits[i][j].split("|")
                qs_ins = qs.split("|")
                # mask low quality insertion as "N"
                cs_replaced = [cs if int(q) >= 10 else "N" for cs, q in zip(cs_ins, qs_ins)]
                # replace to match when if all inserted bases have low quality
                if all(cs == "N" for cs in cs_replaced):
                    cs_replaced = "=" + sequence[j]
                else:
                    cs_replaced = "|".join(cs_replaced)
                cssplits[i][j] = cs_replaced
            elif 0 <= int(qs) < 10:
                cssplits[i][j] = "=" + sequence[j]
    return cssplits


###############################################################################
# Count
###############################################################################


def transpose(cssplit: list[list[str]]):
    return list(zip(*cssplit))


def call_count(cssplit_transposed: list[list[str]]):
    cssplit_counts = []
    for cssplit in cssplit_transposed:
        count = Counter(cssplit)
        count = dict(count)
        cssplit_counts.append(count)
    return cssplit_counts


def call_percent(cssplit_sample_counts):
    cssplit_percent = []
    coverage = sum(cssplit_sample_counts[0].values())
    for counts in cssplit_sample_counts:
        percent = {k: v / coverage * 100 for k, v in counts.items()}
        cssplit_percent.append(percent)
    return cssplit_percent


###############################################################################
# Mask repetitive dels and knock-in
###############################################################################


def pearson_corr(x: list[int], y: list[int]):
    x_diff = np.array(x) - np.mean(x)
    y_diff = np.array(y) - np.mean(y)
    if (np.sqrt(sum(x_diff ** 2)) * np.sqrt(sum(y_diff ** 2))) == 0:
        return 0
    return np.dot(x_diff, y_diff) / (np.sqrt(sum(x_diff ** 2)) * np.sqrt(sum(y_diff ** 2)))


def count_deletons_in_repaet(control_count, sequence):
    list_deletions = []
    count_deletions = []
    repeat_regrex = r"A{3,}|C{3,}|G{3,}|T{3,}|N{3,}"
    repeat_span = list(x.span() for x in re.finditer(repeat_regrex, sequence))
    repeat_idx = 0
    repeat_start = repeat_span[repeat_idx][0]
    repeat_end = repeat_span[repeat_idx][1]
    for i, counts in enumerate(control_count):
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


def find_repetitive_dels(control_count, sample_count, sequence):
    control_repdels = count_deletons_in_repaet(control_count, sequence)
    sample_repdels = count_deletons_in_repaet(sample_count, sequence)
    repeat_regrex = r"A{4,}|C{4,}|G{4,}|T{4,}|N{4,}"
    repeat_span = list(x.span() for x in re.finditer(repeat_regrex, sequence))
    corrs = [pearson_corr(c, s) for c, s in zip(control_repdels, sample_repdels)]
    highcorrs = [list(range(span[0], span[1])) for corr, span in zip(corrs, repeat_span) if corr > 0.95]
    return set(chain.from_iterable(highcorrs))


def recount_repdels_as_match(sample_count, repeat_dels):
    recount = deepcopy(sample_count)
    for i, count in enumerate(recount):
        if i not in repeat_dels:
            continue
        for key, val in count.items():
            if key.startswith("-"):
                key_del = key
                val_del = val
            elif key.startswith("="):
                key_match = key
        try:
            count[key_match] += val_del
        except KeyError:
            pass
        try:
            del count[key_del]
        except KeyError:
            pass
    return recount


def recount_knockin_loci(control_count, knockin_loci, sequence):
    recount = deepcopy(control_count)
    for i, count in enumerate(recount):
        if i not in knockin_loci:
            continue
        key_match = "=" + sequence[i]
        val_match = 0
        for key, val in count.items():
            if key.startswith("-"):
                key_del = key
                val_del = val
            elif key.startswith("="):
                val_match = val
        count.update({key_match: val_match + val_del})
        del count[key_del]
    return recount


###############################################################################
# Discard errors
###############################################################################


def subtract_percentage(control_percent, sample_percent) -> list[dict]:
    sample_subtracted = []
    for cont, samp in zip(control_percent, sample_percent):
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


def lowess_with_confidence_bounds(x, y, eval_x, N=200, conf_interval=0.99, lowess_kw=None):
    """
    Perform Lowess regression and determine a confidence interval by bootstrap resampling
    """
    # Lowess smoothing
    smoothed = sm.nonparametric.lowess(exog=x, endog=y, xvals=eval_x, **lowess_kw)
    # Perform bootstrap resamplings of the data
    # and  evaluate the smoothing at a fixed set of points
    smoothed_values = np.empty((N, len(eval_x)))
    for i in range(N):
        sample = np.random.choice(len(x), len(x), replace=True)
        sampled_x = x[sample]
        sampled_y = y[sample]
        smoothed_values[i] = sm.nonparametric.lowess(exog=sampled_x, endog=sampled_y, xvals=eval_x, **lowess_kw)
    # Get the confidence interval
    sorted_values = np.sort(smoothed_values, axis=0)
    bound = int(N * (1 - conf_interval) / 2)
    bottom = sorted_values[bound - 1]
    top = sorted_values[-bound]
    return smoothed, bottom, top


def discard_common_dels(control_percent, sample_percent):
    per_mismatch = []
    per_del = []
    for cont in control_percent:
        match = 0
        dels = 0
        for key, val in Counter(cont).items():
            if key.startswith("="):
                match += val
            elif key.startswith("-"):
                dels += val
            mismatch = 100 - match
        per_mismatch.append(mismatch)
        per_del.append(dels)
    # Calcurate Lowess
    eval_x = np.linspace(0, 100, 1000)
    smoothed, bottom, top = lowess_with_confidence_bounds(
        np.array(per_mismatch), np.array(per_del), eval_x, N=200, conf_interval=0.99, lowess_kw={"frac": 0.1}
    )
    # Replace common deletion loci
    sample_percent_discarded = deepcopy(sample_percent)
    for i in range(len(sample_percent)):
        per_del = sum(val for key, val in sample_percent[i].items() if key.startswith("-"))
        per_mismatch = 100 - sum(val for key, val in sample_percent[i].items() if key.startswith("="))
        idx = bisect.bisect(eval_x, per_mismatch)
        if per_del - top[idx] > 0.5:
            continue
        for key, val in sample_percent[i].items():
            if key.startswith("-"):
                del sample_percent_discarded[i][key]
    return sample_percent_discarded


###############################################################################
# Annotgate scores
###############################################################################


def allocate_scores(sample_t, sample_discarded):
    scores = []
    for samp, mutation in zip(sample_t, sample_discarded):
        score = []
        if not mutation:
            # # TODO 全インデックスにスコアを当てるとデバッグがし易い。本番では計算量低減のために以下の2行は削除する。
            score = [0 for _ in range(len(samp))]
            scores.append(score)
            continue
        for key, val in mutation.items():
            for s in samp:
                if s == key:
                    score.append(val)
                else:
                    score.append(0)
        scores.append(score)
    return transpose(scores)


def calculate_scores(cssplits_control, cssplits_sample, qscores, sequence, knockin_loci):
    # (1) replace N and Low quality bases
    cssplits_control_replaced = replace_n_to_match(cssplits_control, sequence)
    cssplits_control_replaced = replace_lowquality_to_match(cssplits_control_replaced, qscores, sequence)
    # (2) Formatting...
    control_t = transpose(cssplits_control_replaced)
    control_count = call_count(control_t)
    # (3) Recount knock-in, repeat deletions, common deletions as match
    ## mask knock-in loci as match (matchではなく、捨てたほうが良い？)
    control_recount = recount_knockin_loci(control_count, knockin_loci, sequence)
    ## mask repeat deletions as match (matchではなく、捨てたほうが良い？)
    repeat_dels = find_repetitive_dels(control_recount, sample_count, sequence)
    control_recount = recount_repdels_as_match(control_recount, repeat_dels)
    sample_recount = recount_repdels_as_match(sample_count, repeat_dels)
    control_percent = call_percent(control_recount)
    sample_percent = call_percent(sample_recount)
    sample_percent_discarded = discard_common_dels(control_percent, sample_percent)
    # (4) Discard...
    sample_subtracted = subtract_percentage(control_percent, sample_percent_discarded)
    sample_discarded = discard_common_error(sample_subtracted, 0.5)
    sample_discarded = discard_match(sample_discarded)
    # (5) Todo: Extract DIFFLOCI, DIFFMUTATIONS, DIFFSCORES
    # (6) Annotate scores
    sample_scores = allocate_scores(sample_t, sample_discarded)
    return sample_scores


for i, s in enumerate(sample_discarded):
    if s:
        print(i, s)


def return_labels(scores: list[list[float]], threads: int = 1) -> list[int]:
    labels = OPTICS(min_samples=0.005, n_jobs=threads).fit_predict(scores)
    labels[labels == -1] = np.max(labels) + 1
    return labels.tolist()


# !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

allele = "control"
sv = False
sequence = FASTA_ALLELE[allele]
knockin_loci = KNOCKIN_LOCI[allele]

midsv_control = midsv.read_jsonl((Path(TEMPDIR, "midsv", f"{CONTROL_NAME}_{allele}.jsonl")))
cssplits_control = [cs["CSSPLIT"].split(",") for cs in midsv_control]
qscores = [cs["QSCORE"].split(",") for cs in midsv_control]
cssplits_control = replace_n_to_match(cssplits_control, sequence)
cssplits_control = replace_lowquality_to_match(cssplits_control, qscores, sequence)
