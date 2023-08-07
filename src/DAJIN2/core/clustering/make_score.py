from __future__ import annotations

from typing import Generator
from itertools import groupby
from collections import Counter

from DAJIN2.core.clustering.make_kmer import generate_mutation_kmers


def _call_count(cssplits_sample: Generator[list[str]]) -> list[dict[str, int]]:
    count_kmer = []
    for cssplits_transposed in zip(*cssplits_sample):
        count_kmer.append(dict(Counter(cssplits_transposed)))
    return count_kmer


def _call_percent(counts: list[dict[str, int]]) -> list[dict[str, float]]:
    cssplit_percent = []
    coverage = sum(counts[0].values())
    for count in counts:
        percent = {k: round(v / coverage * 100, 3) for k, v in count.items()}
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
# Handling insertions
###############################################################################


def _group_consecutive_insertions(mutation_loci: dict[str, set[int]]) -> list[tuple[int]]:
    index = sorted(i for i, m in enumerate(mutation_loci) if "+" in m)
    index_grouped = []
    for _, group in groupby(enumerate(index), lambda i_x: i_x[0] - i_x[1]):
        items = [value for _, value in group]
        index_grouped.append(tuple(items))
    return index_grouped


def update_insertion_score(percent_discarded, mutation_loci):
    index_grouped = _group_consecutive_insertions(mutation_loci)
    for index_insertion in index_grouped:
        max_score = 0
        for i in index_insertion:
            for key, val in percent_discarded[i].items():
                if key.split(",")[1].startswith("+"):
                    max_score = max(max_score, val)
        for i in index_insertion:
            for key in percent_discarded[i]:
                if key.split(",")[1].startswith("+"):
                    percent_discarded[i][key] = max_score
    return percent_discarded


###############################################################################
# main
###############################################################################


def make_score(path_sample, path_control, mutation_loci: set[int], knockin_loci: set[int]) -> list[dict[str, float]]:
    cssplits_sample = generate_mutation_kmers(path_sample, mutation_loci, compress_ins=True)
    cssplits_control = generate_mutation_kmers(path_control, mutation_loci, compress_ins=True)
    counts_sample = _call_count(cssplits_sample)
    counts_control = _call_count(cssplits_control)
    percent_sample = _call_percent(counts_sample)
    percent_control = _call_percent(counts_control)
    percent_subtraction = _subtract_percentage(percent_control, percent_sample, knockin_loci)
    percent_discarded = _discard_common_error(percent_subtraction, 0.5)
    percent_discarded = update_insertion_score(percent_discarded, mutation_loci)
    mutation_score = _discard_match_and_n(percent_discarded)
    return mutation_score
