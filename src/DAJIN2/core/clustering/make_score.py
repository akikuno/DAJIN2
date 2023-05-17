from __future__ import annotations

from collections import Counter
from typing import Generator


def _call_count(cssplits_sample: Generator[list[str]], mutation_loci) -> list[dict[str, int]]:
    count_kmer = []
    for i, cssplits_transposed in enumerate(zip(*cssplits_sample)):
        count_kmer.append(dict(Counter(cssplits_transposed)))
    return count_kmer


def _call_percent(counts: list[dict[str, int]]) -> list[dict[str, float]]:
    cssplit_percent = []
    coverage = sum(counts[0].values())
    for count in counts:
        percent = {k: v / coverage * 100 for k, v in count.items()}
        cssplit_percent.append(percent)
    return cssplit_percent


def _subtract_percentage(percent_control, percent_sample, knockin_loci: set(int)) -> list[dict]:
    sample_subtracted = []
    for i, (cont, samp) in enumerate(zip(percent_control, percent_sample)):
        if i in knockin_loci:
            sample_subtracted.append(dict(samp))
            continue
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
    mutation_score_discarded = [dict() for _ in range(len(percent_discarded))]
    for i, mutation_percent in enumerate(percent_discarded):
        if mutation_percent == {} or "" in mutation_percent:
            continue
        for mutation, percent in mutation_percent.items():
            mutation_center = mutation.split(",")[1]
            if mutation_center.startswith("=") or mutation_center == ("N"):
                continue
            mutation_score_discarded[i].update({mutation: percent})
    return mutation_score_discarded


###############################################################################
# main
###############################################################################


def make_score(
    cssplits_sample, cssplits_control, mutation_loci: dict[str, set[int]], knockin_loci: set(int)
) -> list[dict[str, float]]:
    counts_sample = _call_count(cssplits_sample, mutation_loci)
    counts_control = _call_count(cssplits_control, mutation_loci)
    percent_sample = _call_percent(counts_sample)
    percent_control = _call_percent(counts_control)
    percent_subtraction = _subtract_percentage(percent_control, percent_sample, knockin_loci)
    percent_discarded = _discard_common_error(percent_subtraction, 0.5)
    mutation_score = _discard_match_and_n(percent_discarded)
    return mutation_score
