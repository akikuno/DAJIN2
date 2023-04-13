from __future__ import annotations

from collections import Counter


def _call_count(cssplits: list[list[str]]) -> list[dict[str, int]]:
    count_kmer = []
    for cs in list(zip(*cssplits)):
        count_kmer.append(dict(Counter(cs)))
    return count_kmer


def _call_percent(counts: list[dict[str:int]]) -> list[dict[str:float]]:
    cssplit_percent = []
    coverage = sum(counts[0].values())
    for count in counts:
        percent = {k: v / coverage * 100 for k, v in count.items()}
        cssplit_percent.append(percent)
    return cssplit_percent


def _subtract_percentage(percent_control, percent_sample) -> list[dict]:
    sample_subtracted = []
    for cont, samp in zip(percent_control, percent_sample):
        samp = Counter(samp)
        samp.subtract(Counter(cont))
        sample_subtracted.append(dict(samp))
    return sample_subtracted


def _discard_common_error(percent_subtraction, threshold=0.5) -> list[dict]:
    percent_discarded = []
    for samp in percent_subtraction:
        remained = {k: v for k, v in samp.items() if v > threshold}
        percent_discarded.append(remained)
    return percent_discarded


def _discard_match_and_n(percent_discarded) -> list[dict]:
    mutation_score = []
    for samp in percent_discarded:
        if samp == {} or "" in samp:
            mutation_score.append({})
            continue
        cs_center = list(samp.keys())[0].split(",")[1]
        if cs_center.startswith("=") or cs_center == ("N"):
            mutation_score.append({})
            continue
        score = {k: v for k, v in samp.items()}
        mutation_score.append(score)
    return mutation_score


###############################################################################
# main
###############################################################################


def make_score(cssplits_control, cssplits_sample) -> list[dict[str, float]]:
    counts_control = _call_count(cssplits_control)
    counts_sample = _call_count(cssplits_sample)
    percent_control = _call_percent(counts_control)
    percent_sample = _call_percent(counts_sample)
    percent_subtraction = _subtract_percentage(percent_control, percent_sample)
    percent_discarded = _discard_common_error(percent_subtraction, 0.5)
    mutation_score = _discard_match_and_n(percent_discarded)
    return mutation_score
