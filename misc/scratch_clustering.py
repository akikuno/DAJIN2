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
from sklearn.cluster import MeanShift
from scipy.spatial.distance import cosine
from itertools import groupby

# *Stx2 deletion: 2012-2739 (727 bases)
# tests/data/knockout/test_barcode25.fq.gz
# tests/data/knockout/test_barcode30.fq.gz
# tests/data/knockout/design_stx2.fa


def transpose(cssplits):
    return [list(cs) for cs in zip(*cssplits)]


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
            count_del = sum(val for key, val in counts.items() if key.startswith("-")) + 1
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


def sampling(cnt: Counter, size: int):
    elements = []
    probs = []
    coverage = sum(cnt.values())
    for key, val in cnt.items():
        elements.append(key)
        probs.append(val / coverage)
    np.random.seed(1)
    samples = np.random.choice(a=elements, size=size, p=probs)
    return samples


def replace_repdels(transposed_cssplits: list[list[str]], repeat_dels: set):
    replased_cssplits = deepcopy(transposed_cssplits)
    for i, cssplits in enumerate(replased_cssplits):
        if i not in repeat_dels:
            continue
        cnt = Counter(cssplits)
        size = sum(1 for cs in cssplits if cs.startswith("-"))
        if size == 0:
            continue
        for key in list(cnt.keys()):
            if key.startswith("-"):
                del cnt[key]
        samples = sampling(cnt, size)
        iter_samples = iter(samples)
        for j, cs in enumerate(cssplits):
            if cs.startswith("-"):
                replased_cssplits[i][j] = next(iter_samples)
    return replased_cssplits


def replace_both_ends_n(transposed_cssplits: list[list[str]]):
    d_samples = defaultdict(iter)
    for i, cssplits in enumerate(transposed_cssplits):
        cnt = Counter(cssplits)
        if cnt["N"] == 0:  # No N
            continue
        if len(cnt) == 1 and cnt["N"] > 0:  # All N
            continue
        size = cnt["N"]
        del cnt["N"]
        samples = sampling(cnt, size)
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
    return transpose(scores)


###############################################################################
from src.DAJIN2.core import clustering

allele = "control"
sv = True
sequence = DICT_ALLELE[allele]
knockin_loci = KNOCKIN_LOCI[allele]

# Control
midsv_control = midsv.read_jsonl((Path(TEMPDIR, "midsv", f"{CONTROL_NAME}_{allele}.jsonl")))
cssplits_control = [cs["CSSPLIT"].split(",") for cs in midsv_control]
# Sample
cssplits_sample = [cs["CSSPLIT"].split(",") for cs in classif_sample if cs["ALLELE"] == allele and cs["SV"] == sv]
cssplits_corrected_control, cssplits_corrected_sample = clustering.correct_cssplits(
    cssplits_control, cssplits_sample, sequence
)
mutation_score = clustering.make_score(cssplits_corrected_control, cssplits_corrected_sample)
scores_control = clustering.annotate_score(cssplits_corrected_control, mutation_score)
scores_sample = clustering.annotate_score(cssplits_corrected_sample, mutation_score)
scores = scores_sample + scores_control[:1000]
labels = clustering.return_labels(scores, THREADS)
labels_control = labels[len(scores_sample) :]
labels_sample = labels[: len(scores_sample)]
labels_merged = clustering.merge_clusters(labels_control, labels_sample)


def add_labels(classif_sample, CONTROL_NAME, DICT_ALLELE, KNOCKIN_LOCI, TEMPDIR, THREADS: int):
    labels = []
    max_label = 0
    classif_sample.sort(key=lambda x: (x["ALLELE"], x["SV"]))
    for (allele, sv), group in groupby(classif_sample, key=lambda x: (x["ALLELE"], x["SV"])):
        sequence = DICT_ALLELE[allele]
        knockin_loci = KNOCKIN_LOCI[allele]
        # Control
        midsv_control = midsv.read_jsonl((Path(TEMPDIR, "midsv", f"{CONTROL_NAME}_{allele}.jsonl")))
        cssplits_control = [cs["CSSPLIT"].split(",") for cs in midsv_control]
        # Sample
        cssplits_sample = [cs["CSSPLIT"].split(",") for cs in group if cs["ALLELE"] == allele and cs["SV"] == sv]
        cssplits_corrected_control, cssplits_corrected_sample = clustering.correct_cssplits(
            cssplits_control, cssplits_sample, sequence
        )
        mutation_score = clustering.make_score(cssplits_corrected_control, cssplits_corrected_sample)
        scores_control = clustering.annotate_score(cssplits_corrected_control, mutation_score)
        scores_sample = clustering.annotate_score(cssplits_corrected_sample, mutation_score)
        scores = scores_sample + scores_control[:1000]
        labels = clustering.return_labels(scores, THREADS)
        labels_control = labels[len(scores_sample) :]
        labels_sample = labels[: len(scores_sample)]
        labels_merged = clustering.merge_clusters(labels_control, labels_sample)
        labels_reorder = clustering.reorder_labels(labels_merged, start=max_label)
        max_label = max(max_label, labels_reorder)
        labels.append(labels_reorder)
    clust_sample = deepcopy(classif_sample)
    for clust, label in zip(clust_sample, labels):
        clust["LABEL"] = label
    return clust_sample


# Mask common mutations
# transpose_control = transpose(cssplits_corrected_control)
# transpose_sample = transpose(cssplits_corrected_sample)
# count_compensated_control = call_count(transpose_control)
# count_compensated_sample = call_count(transpose_sample)
# percent_control = call_percent(count_compensated_control)
# percent_sample = call_percent(count_compensated_sample)
# percent_subtraction = subtract_percentage(percent_control, percent_sample)
# percent_discarded = discard_common_error(percent_subtraction, 0.5)
# mutation_score = discard_match(percent_discarded)

# scores_control = annotate_scores(transpose_control, mutation_score)
# scores_sample = annotate_scores(transpose_sample, mutation_score)

# Clustering
