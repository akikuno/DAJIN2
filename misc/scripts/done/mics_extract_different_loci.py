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


np.random.seed(1)

###############################################################################
# Replace N and Low Quality
###############################################################################


def replace_n_to_match(cssplits: list[list[str]], sequence: str) -> list[str]:
    cssplits_replace = deepcopy(cssplits)
    for idx, cssplit in enumerate(cssplits_replace):
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
        cssplits_replace[idx] = cssplit[::-1]
    return cssplits_replace


def replace_lowquality_to_match(cssplits: list[list[str]], qscores: list[list[int]], sequence: str) -> list[str]:
    cssplits_replace = deepcopy(cssplits)
    for i, qscore in enumerate(qscores):
        for j, qs in enumerate(qscore):
            if "|" in qs:
                cs_ins = cssplits_replace[i][j].split("|")
                qs_ins = qs.split("|")
                # mask low quality insertion as "N"
                cs_replaced = [cs if int(q) >= 10 else "N" for cs, q in zip(cs_ins, qs_ins)]
                # replace to match when if all inserted bases have low quality
                if all(cs == "N" for cs in cs_replaced):
                    cs_replaced = "=" + sequence[j]
                else:
                    cs_replaced = "|".join(cs_replaced)
                cssplits_replace[i][j] = cs_replaced
            elif 0 <= int(qs) < 10:
                cssplits_replace[i][j] = "=" + sequence[j]
    return cssplits_replace


###############################################################################
# Count
###############################################################################


def transpose(cssplit: list[list[str]]) -> list[list[str]]:
    return list(zip(*cssplit))


def call_count(cssplit_transposed: list[list[str]]) -> list[dict[str:int]]:
    cssplit_counts = []
    for cssplit in cssplit_transposed:
        count = Counter(cssplit)
        count = dict(count)
        cssplit_counts.append(count)
    return cssplit_counts


def call_percent(cssplit_counts: list[dict[str:int]]) -> list[dict[str:int]]:
    cssplit_percent = []
    coverage = sum(cssplit_counts[0].values())
    for counts in cssplit_counts:
        percent = {k: v / coverage * 100 for k, v in counts.items()}
        cssplit_percent.append(percent)
    return cssplit_percent


###############################################################################
# Mask repetitive dels and knock-in
###############################################################################


def pearson_corr(x: list[int], y: list[int]) -> int:
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


def find_repetitive_dels(control_count, sample_count, sequence) -> set[int]:
    control_repdels = count_deletons_in_repaet(control_count, sequence)
    sample_repdels = count_deletons_in_repaet(sample_count, sequence)
    repeat_regrex = r"A{3,}|C{3,}|G{3,}|T{3,}|N{3,}"
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


# ------------------------------------------------------------------------------
# Find knock-in loci
# ------------------------------------------------------------------------------


def find_knockin_loci(cssplits_count, sequence) -> set[int]:
    """
    Alignments between control and knock-in alleles produses deletion loci in control.
    The deletion loci are knock-in loci, not sequencing error, so they should be ignored.
    """
    if len(cssplits_count) <= len(sequence):
        return set()
    coverage = sum(cssplits_count[0].values())
    cssplit_mostcommon = [Counter(cs).most_common()[0] for cs in cssplits_count]
    knockin_loci = set()
    for i, (key, val) in enumerate(cssplit_mostcommon):
        val_percent = val / coverage * 100
        if key.startswith("-") and val_percent > 80:
            knockin_loci.add(i)
    return knockin_loci


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
    return transpose(scores)


###############################################################################
# Summarize
###############################################################################

# def return_labels(scores: list[list[float]], threads: int = 1) -> list[int]:
#     labels = OPTICS(min_samples=0.005, n_jobs=threads).fit_predict(scores)
#     labels[labels == -1] = np.max(labels) + 1
#     return labels.tolist()


# def call_consensus(cssplits_control, sequence):
#     cssplits_transpose = transpose(cssplits_control)
#     repeat_regrex = r"A{3,}|C{3,}|G{3,}|T{3,}|N{3,}"
#     repeat_span = list(x.span() for x in re.finditer(repeat_regrex, sequence))
#     repeat_loci = set()
#     for s, e in repeat_span:
#         for i in list(range(s, e)):
#             repeat_loci.add(i)
#     consensus = []
#     for i, cssplits in enumerate(cssplits_transpose):
#         percent = defaultdict(int)
#         coverage = len(cssplits) - Counter(cssplits)["N"]
#         for cs in cssplits:
#             if cs.startswith("+"):
#                 percent["+"] += 1 / coverage * 100
#             else:
#                 percent[cs] += 1 / coverage * 100
#         cons, per = sorted(percent.items(), key=lambda x: -x[1])[0]
#         if i in repeat_loci and cons.startswith("-") and not per > 90:
#             cons = "=" + sequence[i]
#         elif per < 50:
#             cons = "=" + sequence[i]
#         consensus.append(cons)
#     return consensus


# x = call_consensus(cssplits_control, sequence)
# [(i, xx) for i, xx in enumerate(x) if not xx.startswith("=")]
# *=============================================================================================

allele = "deletion"
sv = True
sequence = FASTA_ALLELE[allele]
knockin_loci = KNOCKIN_LOCI[allele]

# Control
midsv_control = midsv.read_jsonl((Path(TEMPDIR, "midsv", f"{CONTROL_NAME}_{allele}.jsonl")))
cssplits_control = [cs["CSSPLIT"].split(",") for cs in midsv_control]
transpose_control = transpose(cssplits_control)
count_control = call_count(transpose_control)


qscores = [cs["QSCORE"].split(",") for cs in midsv_control]

cssplits_replaced = replace_n_to_match(cssplits_control, sequence)
cssplits_replaced = replace_lowquality_to_match(cssplits_replaced, qscores, sequence)
transpose_control = transpose(cssplits_replaced)
count_control = call_count(transpose_control)

knockin_loci = find_knockin_loci(count_control, sequence)
count_control = recount_knockin_loci(count_control, knockin_loci, sequence)

# Sample
cssplits_sample = [cs["CSSPLIT"].split(",") for cs in classif_sample if cs["ALLELE"] == allele and cs["SV"] == sv]
transpose_sample = transpose(cssplits_sample)
qscores = [cs["QSCORE"].split(",") for cs in classif_sample if cs["ALLELE"] == allele and cs["SV"] == sv]
cssplits_replaced = replace_n_to_match(cssplits_sample, sequence)
cssplits_replaced = replace_lowquality_to_match(cssplits_replaced, qscores, sequence)
transpose_sample = transpose(cssplits_replaced)
count_sample = call_count(transpose_sample)

[cs for cs in cssplits_sample[-9] if re.search(r"[acgtn]", cs)]
cs = ",".join(cssplits_sample[-9])
for classif in classif_sample:
    if classif["CSSPLIT"] == cs:
        print(classif["QNAME"])

transpose_sample[2100][-9]
# ControlもSampleも使用している
repeat_dels = find_repetitive_dels(count_control, count_sample, sequence)  # *依存

recount_control = recount_repdels_as_match(count_control, repeat_dels)
recount_sample = recount_repdels_as_match(count_sample, repeat_dels)
percent_control = call_percent(recount_control)
percent_sample = call_percent(recount_sample)
percent_discard_sample = discard_common_dels(percent_control, percent_sample)  # *依存
# (4) Discard...
percent_discard_control = discard_common_error(percent_control, 0.5)
percent_discard_control = discard_match(percent_discard_control)
percent_subtract_sample = subtract_percentage(percent_control, percent_discard_sample)
percent_discard_sample = discard_common_error(percent_subtract_sample, 0.5)
percent_discard_sample = discard_match(percent_discard_sample)

# (5) Todo: Extract DIFFLOCI, DIFFMUTATIONS, DIFFSCORES
# (6) Annotate scores
scores_control = allocate_scores(transpose_control, percent_discard_sample)
scores_sample = allocate_scores(transpose_sample, percent_discard_sample)

scores_control_subset = scores_control[:1000]
scores = scores_sample + scores_control_subset
labels = MeanShift(n_jobs=THREADS).fit(scores).labels_
Counter(labels)

labels_control = deepcopy(labels)[len(scores_sample) :].tolist()
labels_sample = deepcopy(labels)[: len(scores_sample)].tolist()


def merge_mixed_cluster(labels_control, labels_sample):
    labels_all = labels_sample + labels_control
    labels_merged = deepcopy(labels_sample)
    coverage_control = len(labels_control)
    labels_percent_control = [0] * len(np.unique(labels_all))
    for i, label in enumerate(labels_control):
        labels_percent_control[label] += 1 / coverage_control * 100
    labels_mixed = {i for i, per in enumerate(labels_percent_control) if per > 0.5}
    max_label = max(labels_merged) + 1
    for i, label in enumerate(labels_sample):
        if label in labels_mixed:
            labels_merged[i] = max_label
    return labels_merged


labels_sample_merge = {
    i for i, v in Counter(labels_sample).items() if v / sum(Counter(labels_sample).values()) * 100 < 0.5
}
max_label = max(labels_sample) + 1
for i, label in enumerate(labels_sample):
    if label in labels_sample_merge:
        labels_sample[i] = max_label

Counter(labels_sample)


def extract_apparent_mutations(cssplits_sample, labels_sample) -> defaultdict[list[str]]:
    labels_diffloci = defaultdict(list)
    for label_num in Counter(labels_sample).keys():
        cssplits_label = [score for score, label in zip(cssplits_sample, labels_sample) if label == label_num]
        readnum = len(cssplits_label)
        transposed = list(zip(*cssplits_label))
        positions = []
        for i, (s, p) in enumerate(zip(transposed, percent_discarded)):
            if not p:
                continue
            if Counter(s).most_common()[0][1] / readnum * 100 > 80:
                key = Counter(s).most_common()[0][0]
                if not key.startswith("="):
                    positions.append(i)
        labels_diffloci[label_num] = " ".join(map(str, positions))
    return labels_diffloci


labels_diffloci = extract_diffloci(cssplits_sample, labels_sample)
labels_diffloci_unique = set(labels_diffloci.values())

d = {label: i + 1 for i, label in enumerate(labels_diffloci_unique)}
x = deepcopy(labels_diffloci)
for label, loci in labels_diffloci.items():
    x[label] = d[loci]

for i, label in enumerate(labels_sample):
    labels_sample[i] = x[label]

Counter(labels_sample)


# Reorder
def order_labels(labels: list[int]) -> list[int]:
    labels_ordered = deepcopy(labels)
    num = 0
    d = defaultdict(int)
    for i, l in enumerate(labels_ordered):
        if not d[l]:
            num += 1
            d[l] = num
        labels_ordered[i] = d[l]
    return labels_ordered


Counter(labels_sample)


# d = defaultdict(int)
# d_all = defaultdict(int)
# for cs, label in zip(cssplits_sample, labels_sample):
#     if label == 2:
#         d[cs[828]] += 1
#     d_all[cs[828]] += 1

# d_all
# count = 0
# count_non = 0
# for s in scores_label:
#     if s[828] > 0:
#         count += 1
#     else:
#         count_non += 1

# count
# count_non

# cons_results = []
# for label_num in sorted(list(Counter(labels_sample).keys())):
#     cons = dict()
#     cons.update({"LABEL": label_num})
#     cssplits_label = [score for score, label in zip(cssplits_sample, labels_sample) if label == label_num]
#     readnum = len(cssplits_label)
#     cons.update({"READNUM": readnum})
#     transposed = transpose(cssplits_label)
#     positions = []
#     for i, (s, p) in enumerate(zip(transposed, percent_discard_sample)):
#         if not p:
#             continue
#         if len(Counter(s)) == 1 or Counter(s).most_common()[0][1] / readnum * 100 > 80:
#             key = Counter(s).most_common()[0][0]
#             if key.startswith("="):
#                 continue
#         else:
#             continue
#         positions.append([i, key])
#     cons.update({"LOCI": positions})
#     cons_results.append(cons)

# cons_results.sort(key=lambda x: x["READNUM"])
# pprint(cons_results)

# i = 2755
# sample_percent[i]
# control_percent[i]


# label_num = 18
# x = [cs for cs, label in zip(cssplits_sample, labels) if label == label_num]
# transpose(x)[2475]
# sample_percent[2475]
# control_percent[2475]
# Counter(transpose(cssplits_control)[2475])
# Counter(transpose(cssplits_sample)[2475])
# from collections import defaultdict

# def call_consensus(cssplits_sample_zero, sequence):
#     cssplits_t = transpose(cssplits_sample_zero)
#     coverage = len(cssplits_sample_zero)
#     consensus = []
#     for i, cssplits in enumerate(cssplits_t):
#         percent = defaultdict(int)
#         for cs in cssplits:
#             percent[cs] += 1 / coverage * 100
#         cons, per = sorted(percent.items(), key=lambda x: -x[1])[0]
#         if per < 80:
#             cons = "=" + sequence[i]
#         consensus.append(cons)
#     return consensus


# for i, cs in enumerate(call_consensus(cssplits_sample_zero, sequence)):
#     if not cs.startswith("="):
#         print(i, cs)


# # ? >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# """
# コントロールとサンプルの変異率をペアで比べることで、エラーと真の変異の区別をつける？
# """

# control_percent[2721]
# sample_percent[2721]

# # def recount_control_error(control_percent, sample_percent, threshold=1):
# paired_percent = []
# for cont, samp in zip(control_percent, sample_percent):
#     keys = cont.keys() | samp.keys()
#     for key in keys:
#         if key.startswith("="):
#             continue
#         # if not (4 < val_cont < 5):
#         #     continue
#         try:
#             cont_val = cont[key]
#         except KeyError:
#             cont.update({key: 0})
#         try:
#             samp_val = samp[key]
#         except KeyError:
#             continue
#         paired_percent.append((key, cont_val, samp_val))

# # for cont, samp in zip(control_percent, sample_percent):
# #     for key_cont, val_cont in cont.items():
# #         if key_cont.startswith("="):
# #             continue
# #         # if not (4 < val_cont < 5):
# #         #     continue
# #         try:
# #             val_samp = samp[key_cont]
# #         except KeyError:
# #             continue
# #         paired_percent.append((key_cont, val_cont, val_samp))


# paired_percent.sort(key=lambda x: [x[1], -x[2]])
# paired_percent[:10]

# x = [x for _, x, y in paired_percent]
# y = [y for _, x, y in paired_percent]

# import matplotlib.pyplot as plt
# import seaborn as sns

# # sns.regplot(x, y, color="blue")
# plt.plot(x, y, "o")
# plt.show()

# # # X = [[x, y] for _, x, y in paired_percent]

# # from sklearn.svm import OneClassSVM

# # clf = OneClassSVM(gamma="auto").fit(X)
# # labels = clf.predict(X)
# # from collections import defaultdict

# # d = defaultdict(int)
# # for l in labels:
# #     d[l] += 1

# # d
# # len(x)
# # x[:10]
# # y[:10]
# # import matplotlib.pyplot as plt

# # plt.plot(x, y, "o")
# # plt.show()

# for i, samp in enumerate(sample_percent):
#     for val_cont in samp.values():
#         if val_cont == 16.10486891385768:
#             print(i)

# # i = 1781
# # pprint(control_percent[i - 5 : i + 5])
# # sample_percent[i - 5 : i + 5]

# # i = 1000
# # cont, samp = control_percent[i], sample_percent[i]

# # ? <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
