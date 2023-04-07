from __future__ import annotations

from collections import Counter, defaultdict


def call_count(cssplits: list[list[str]]) -> list[dict[str, int]]:
    """Count cssplits within 3-mer range.
    Args:
        cssplits (list[list[str]])
    Returns:
        list[dict[str, int]]: Both ends are counted as "N" to keep sequence length.
    """
    count_kmer = defaultdict(Counter)
    for cssplit in cssplits:
        for i in range(1, len(cssplit) - 1):
            kmer = ",".join([cssplit[i - 1], cssplit[i], cssplit[i + 1]])
            count_kmer[i] += Counter([kmer])
    coverage = len(cssplits)
    count_score = [{"N,N,N": coverage}]
    count_score += [dict(count_kmer[i]) for i in range(1, len(cssplit) - 1)]
    count_score += [{"N,N,N": coverage}]
    return count_score


def call_percent(counts: list[dict[str:int]]) -> list[dict[str:float]]:
    cssplit_percent = []
    coverage = sum(counts[0].values())
    for count in counts:
        percent = {k: v / coverage * 100 for k, v in count.items()}
        cssplit_percent.append(percent)
    return cssplit_percent


def subtract_percentage(percent_control, percent_sample, knockin_loci) -> list[dict]:
    sample_subtracted = []
    for i, (cont, samp) in enumerate(zip(percent_control, percent_sample)):
        if i in knockin_loci:
            sample_subtracted.append(samp)
            continue
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
        remained = {k: v for k, v in samp.items() if not k.split(",")[1].startswith("=")}
        sample_discarded.append(remained)
    return sample_discarded


###############################################################################
# main
###############################################################################


def make_score(cssplits_control, cssplits_sample, knockin_loci):
    counts_control = call_count(cssplits_control)
    counts_sample = call_count(cssplits_sample)
    percent_control = call_percent(counts_control)
    percent_sample = call_percent(counts_sample)
    percent_subtraction = subtract_percentage(percent_control, percent_sample, knockin_loci)
    percent_discarded = discard_common_error(percent_subtraction, 0.5)
    mutation_score = discard_match(percent_discarded)
    return mutation_score
