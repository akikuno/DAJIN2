from __future__ import annotations

from collections import Counter
from collections.abc import Iterator
from itertools import groupby

from DAJIN2.core.clustering.kmer_generator import generate_mutation_kmers


def subset_scores(labels: list[int], scores: list[int], label_most: int, length: int = 1000) -> list[int]:
    subset = [score for label, score in zip(labels, scores) if label == label_most]
    return subset[:length]


def call_count(cssplits: Iterator[list[str]]) -> list[dict[str, int]]:
    return [dict(Counter(cssplit)) for cssplit in zip(*cssplits)]


def call_percent(counts: list[dict[str, int]]) -> list[dict[str, float]]:
    coverage = sum(counts[0].values())
    return [{k: round(v / coverage * 100, 3) for k, v in count.items()} for count in counts]


def subtract_percentage(
    percent_sample: list[dict[str, float]], percent_control: list[dict[str, float]], knockin_loci: set[int]
) -> list[dict[str, float]]:
    """Note: The subtraction operation with collections.Counter class removes elements from the result that have a value of 0 or negative."""
    return [
        dict(samp) if i in knockin_loci else dict(Counter(samp) - Counter(cont))
        for i, (samp, cont) in enumerate(zip(percent_sample, percent_control))
    ]


def discard_common_error(percent_subtraction: list[dict[str, float]], threshold: float = 0.5) -> list[dict]:
    return [{k: v for k, v in samp.items() if v > threshold} for samp in percent_subtraction]


def discard_matches_and_ns(percent_discarded: list[dict[str, float]]) -> list[dict[str, float]]:
    result = []
    for mutation_percent in percent_discarded:
        filtered_dict = {}
        for mutation, percent in mutation_percent.items():
            # Extract the center part of the mutation string
            mutation_center = mutation.split(",")[1]
            # Skip mutations that start with '=' or are 'N'
            if mutation_center.startswith("=") or mutation_center == "N":
                continue
            # Add to the filtered dictionary if it passed the conditions
            filtered_dict[mutation] = percent
        result.append(filtered_dict)

    return result


###############################################################################
# Handling insertions
###############################################################################


def group_consecutive_insertions(mutation_loci: dict[str, set[int]]) -> list[tuple[int]]:
    index = sorted(i for i, m in enumerate(mutation_loci) if "+" in m)
    index_grouped = []
    for _, group in groupby(enumerate(index), lambda i_x: i_x[0] - i_x[1]):
        items = [value for _, value in group]
        index_grouped.append(tuple(items))
    return index_grouped


def find_max_insertion_score_in_group(percent_discarded: list[dict[str, float]], index_group: tuple[int]) -> float:
    max_score = 0
    for i in index_group:
        for key, val in percent_discarded[i].items():
            if key.split(",")[1].startswith("+"):
                max_score = max(max_score, val)
    return max_score


def update_insertion_scores_in_group_to_max(
    percent_discarded: list[dict[str, float]], index_group: tuple[int], max_score: float
) -> list[dict[str, float]]:
    for i in index_group:
        for key in percent_discarded[i]:
            if key.split(",")[1].startswith("+"):
                percent_discarded[i][key] = max_score
    return percent_discarded


def update_insertion_score(
    percent_discarded: list[dict[str, float]], mutation_loci: dict[str, set[int]]
) -> list[dict[str, float]]:
    """Update the insertion scores in the percent_discarded list based on the max score in each group of consecutive insertions.

    This function groups consecutive insertion indices and updates the scores of insertions within those groups to the maximum score found in each group.

    Args:
        percent_discarded (list[dict[str, float]]): A list of dictionaries. Each dictionary contains the scores of mutations at a particular locus, represented as a percentage.
        mutation_loci (dict[str, set[int]]): A dictionary mapping mutation types to sets of indices where those mutations occur.

    Returns:
        list[dict[str, float]]: The updated list of dictionaries containing the new insertion scores.

    """
    index_grouped = group_consecutive_insertions(mutation_loci)
    for index_insertion in index_grouped:
        max_score = find_max_insertion_score_in_group(percent_discarded, index_insertion)
        percent_discarded = update_insertion_scores_in_group_to_max(percent_discarded, index_insertion, max_score)
    return percent_discarded


###############################################################################
# make_score
###############################################################################


def make_score(
    path_sample, path_control, mutation_loci: list[set[int]], knockin_loci: list[set[int]]
) -> list[dict[str, float]]:
    cssplits_sample = generate_mutation_kmers(path_sample, mutation_loci, compress_ins=True)
    cssplits_control = generate_mutation_kmers(path_control, mutation_loci, compress_ins=True)
    counts_sample = call_count(cssplits_sample)
    counts_control = call_count(cssplits_control)
    percent_sample = call_percent(counts_sample)
    percent_control = call_percent(counts_control)
    percent_subtraction = subtract_percentage(percent_sample, percent_control, knockin_loci)
    percent_discarded = discard_common_error(percent_subtraction, 0.5)
    percent_discarded = update_insertion_score(percent_discarded, mutation_loci)
    mutation_score = discard_matches_and_ns(percent_discarded)
    return mutation_score


###############################################################################
# annotate_score
###############################################################################


def annotate_score(path_sample, mutation_score, mutation_loci, is_control=False) -> Iterator[list[float]]:
    for cssplit_kmer in generate_mutation_kmers(path_sample, mutation_loci):
        score = [0 for _ in range(len(cssplit_kmer))]
        for i, (cs_kmer, mut_score) in enumerate(zip(cssplit_kmer, mutation_score)):
            if mut_score == {}:
                continue
            # Mutation sites are not considered in controls because they should be sample-specific.
            if is_control and cs_kmer.split(",")[1][0] in mutation_loci[i]:
                continue
            if cs_kmer in mut_score:
                score[i] = mut_score[cs_kmer]
        yield score
