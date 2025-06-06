from __future__ import annotations

from collections import defaultdict
from pathlib import Path

from scipy.stats import fisher_exact

from DAJIN2.utils import io


def convert_nested_defaultdict_to_dict(d: dict | defaultdict) -> dict:
    """Convert nested defaultdict to a regular dict."""
    if isinstance(d, (defaultdict, dict)):
        return {k: convert_nested_defaultdict_to_dict(v) for k, v in d.items()}
    return d


def count_strandless_of_indels(
    path_midsv_sample: Path, anomal_loci_transposed: list[set[str]]
) -> dict[int, dict[str, dict[str, int]]]:
    count_strandness = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))

    for samp in io.read_jsonl(path_midsv_sample):
        for i, cs in enumerate(samp["CSSPLIT"].split(",")):
            if anomal_loci_transposed[i] == set():
                continue
            for mut in anomal_loci_transposed[i]:
                if not cs.startswith(mut):
                    continue
                count_strandness[i][mut][samp["STRAND"]] += 1
                count_strandness[i][mut].setdefault("+", 0)
                count_strandness[i][mut].setdefault("-", 0)
    return convert_nested_defaultdict_to_dict(count_strandness)


def extract_sequence_errors_in_strand_biased_loci(
    path_midsv_sample: Path, anomal_loci_transposed: list[set[str]]
) -> dict[str, set[int]]:
    coverage = io.count_newlines(path_midsv_sample)
    num_strand_plus = sum(True for m in io.read_jsonl(path_midsv_sample) if m["STRAND"] == "+")
    num_strand_minus = coverage - num_strand_plus

    count_strandness = count_strandless_of_indels(path_midsv_sample, anomal_loci_transposed)

    # A sample is considered "biased" if it exhibits a strong positive or negative strandness bias
    # and the bias is evident when compared to the overall sample set.
    bias_in_strandness = {"+": set(), "-": set(), "*": set()}
    for i, d in count_strandness.items():
        for mut, strandness in d.items():
            # numerical bias detection
            plus_ratio = strandness["+"] / (strandness["+"] + strandness["-"])
            is_biased_ratio = False if 0.2 < plus_ratio < 0.8 else True
            # statistical bias detection
            table = [[strandness["+"], strandness["-"]], [num_strand_plus, num_strand_minus]]
            res = fisher_exact(table, alternative="two-sided")

            if is_biased_ratio and res.pvalue < 0.01:
                bias_in_strandness[mut].add(i)

    return bias_in_strandness
