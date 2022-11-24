from __future__ import annotations
from collections import Counter


def transpose(cssplits):
    return [list(cs) for cs in zip(*cssplits)]


def call_count(transpose_cssplits: list[list[str]]) -> list[dict[str:int]]:
    cssplit_counts = []
    for cssplit in transpose_cssplits:
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
# main
###############################################################################


def make_score(cssplits_control, cssplits_sample):
    transpose_cssplits_control = transpose(cssplits_control)
    transpose_cssplits_sample = transpose(cssplits_sample)
    count_control = call_count(transpose_cssplits_control)
    count_sample = call_count(transpose_cssplits_sample)
    percent_control = call_percent(count_control)
    percent_sample = call_percent(count_sample)
    percent_subtraction = subtract_percentage(percent_control, percent_sample)
    percent_discarded = discard_common_error(percent_subtraction, 0.5)
    mutation_score = discard_match(percent_discarded)
    return mutation_score

