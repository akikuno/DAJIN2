from __future__ import annotations

import re
import bisect
from pathlib import Path
from typing import Generator
from collections import defaultdict

from DAJIN2.utils import config

"""
To suppress the following warnings from `scipy.wilcoxon`:
UserWarning: Exact p-value calculation does not work if there are zeros.
"""
config.set_warnings_ignore()


import numpy as np
from sklearn.cluster import MiniBatchKMeans

from DAJIN2.utils import io
from DAJIN2.core.preprocess import homopolymer_handler


def count_indels(midsv_sample: Generator[dict], sequence: str) -> dict[str, list[int]]:
    len_sequence = len(sequence)
    count = {"=": [0] * len_sequence, "+": [0] * len_sequence, "-": [0] * len_sequence, "*": [0] * len_sequence}
    for samp in midsv_sample:
        for i, cs in enumerate(samp["CSSPLIT"].split(",")):
            if cs == "N" or re.search(r"a|c|g|t|n", cs):
                continue
            if cs.startswith("="):
                count["="][i] += 1
            elif cs.startswith("+"):
                count["+"][i] += 1
            elif cs.startswith("-"):
                count["-"][i] += 1
            elif cs.startswith("*"):
                count["*"][i] += 1
    return count


def normalize_indels(count: dict[str, list[int]]) -> dict[str, np.array]:
    count_normalized = dict()
    match_count = np.array(count["="])
    for mut, indel_count in count.items():
        numerator = np.array(indel_count)
        denominator = numerator + match_count
        count_normalized[mut] = np.where(denominator == 0, 0, numerator / denominator)
    return count_normalized


def minimize_mutation_counts(
    indels_control: dict[str, np.array], indels_sample: dict[str, np.array]
) -> dict[str, np.array]:
    """
    In cases where control has a larger value than sample, adjust the value of sample to match that of control.
    """
    indels_control_minimized = dict()
    for mut in {"+", "-", "*"}:
        indels_control_minimized[mut] = np.minimum(indels_control[mut], indels_sample[mut])
    return indels_control_minimized


def split_kmer(indels: dict[str, np.array], kmer: int = 11) -> dict[str, np.array]:
    results = defaultdict(list)
    center = kmer // 2
    for mut, value in indels.items():
        for i in range(len(value)):
            if center <= i <= len(value) - center:
                start = i - center
                if kmer % 2 == 0:
                    end = i + center
                else:
                    end = i + center + 1
                results[mut].append(value[start:end])
            else:
                results[mut].append(np.array([0] * kmer))
    return results


###########################################################
# Extract Anomalous Loci
###########################################################


def detect_anomalies(values_subtract: np.array) -> set[int]:
    """
    Detect anomalies and return indices of outliers.

    It classifies data points into two classes: 1 (inliers) and -1 (outliers).
    However, depending on how the "normal" class was learned, either of these classes
    might represent the true anomalies in the context of this problem.

    This function returns the indices of the class with the higher mean of values_subtract
    values, as this class is considered to be the true anomalies.
    """
    values_subtract_reshaped = values_subtract.reshape(-1, 1)
    kmeans = MiniBatchKMeans(n_clusters=2, random_state=0)
    _ = kmeans.fit_predict(values_subtract_reshaped)
    threshold = kmeans.cluster_centers_.mean()
    return {i for i, v in enumerate(values_subtract_reshaped) if v > threshold}


def extract_anomal_loci(
    indels_normalized_sample,
    indels_normalized_control,
    thresholds: dict[str, float],
) -> dict[str, set[int]]:
    results = dict()
    for mut in {"+", "-", "*"}:
        values_sample = indels_normalized_sample[mut]
        values_control = indels_normalized_control[mut]
        values_subtract = values_sample - values_control
        """"
        When the result of subtraction is threshold (%) or less, ignore it as 0
        """
        threshold = thresholds[mut]
        values_subtract = np.where(values_subtract <= threshold, 0, values_subtract)
        idx_outliers = detect_anomalies(values_subtract)
        results[mut] = idx_outliers
    return results


###########################################################
# Homolopolymer region
###########################################################


def discard_errors_in_homopolymer(loci: dict[str, set[int]], errors: dict[str, set[int]]) -> dict[str, set[int]]:
    """Remove detected errors in homopolymer regions from the candidate loci."""
    return {mut: loci[mut] - errors[mut] for mut in {"+", "-", "*"}}


###########################################################
# Merge contiguous insertions/deletions
###########################################################


def count_elements_within_range(arr, lower_bound, upper_bound):
    """
    Counts the number of elements within a specified range in a sorted array.
    """
    start_index = bisect.bisect_left(arr, lower_bound)
    end_index = bisect.bisect_right(arr, upper_bound)
    return end_index - start_index


def merge_index_of_consecutive_indel(mutation_loci: dict[str, set[int]]) -> dict[str, set[int]]:
    """Treat as contiguous indels if there are insertions/deletions within five bases of each other"""
    mutation_loci_merged = dict()

    """Reflect point mutations as they are"""
    mutation_loci_merged["*"] = mutation_loci["*"]

    """Merge if indels are within 10 bases"""
    for mut in ["+", "-"]:
        idx_indel = sorted(mutation_loci[mut])
        idx_indel_merged = set(idx_indel)
        for i in range(len(idx_indel) - 1):
            idx_1 = idx_indel[i]
            idx_2 = idx_indel[i + 1]
            """If everything from idx_1 to idx_2 is already considered as indels, then skip it."""
            if count_elements_within_range(idx_indel, idx_1 + 1, idx_2 - 1) == idx_2 - idx_1 + 1:
                continue
            if idx_1 + 10 > idx_2:
                for i in range(idx_1 + 1, idx_2):
                    idx_indel_merged.add(i)
        mutation_loci_merged[mut] = idx_indel_merged

    """Additional logic for mutation enrichment within 10 bases on both ends"""
    for mut in ["+", "-"]:
        idx_indel = sorted(mutation_loci_merged[mut])
        idx_indel_merged = set(idx_indel)
        for i in range(len(idx_indel) - 1):
            idx_1 = idx_indel[i]
            idx_2 = idx_indel[i + 1]
            """If everything from idx_1 to idx_2 is already considered as indels, then skip it."""
            if count_elements_within_range(idx_indel, idx_1 + 1, idx_2 - 1) == idx_2 - idx_1 + 1:
                continue
            """If the distance between idx_1 and idx_2 is more than 20 bases, then skip it."""
            if idx_1 + 20 < idx_2:
                continue
            count_left = count_elements_within_range(idx_indel, idx_1 - 11, idx_1 - 1)
            count_right = count_elements_within_range(idx_indel, idx_2 + 1, idx_2 + 11)
            """
            If 8 out of the 10 bases at both ends are indels,
            then everything from idx_1 to idx_2 will be considered as indels.
            """
            if count_left >= 8 and count_right >= 8:
                for i in range(idx_1 + 1, idx_2):
                    idx_indel_merged.add(i)
        mutation_loci_merged[mut] = idx_indel_merged

    return mutation_loci_merged


###########################################################
# main
###########################################################


def summarize_indels(path_midsv: Path, sequence: str) -> tuple:
    """Returns indels, coverages, normalized indels, and kmer indels."""
    indels_count = count_indels(io.read_jsonl(path_midsv), sequence)
    indels_normalized = normalize_indels(indels_count)

    return indels_count, indels_normalized


def merge_loci(dissimilar_loci: dict[str, set], anomal_loci: dict[str, set]) -> dict[str, set]:
    mutation_loci = dict()
    for mut in {"+", "-", "*"}:
        mutation_loci[mut] = dissimilar_loci[mut] | anomal_loci[mut]
    return mutation_loci


def add_knockin_loci(candidate_loci: dict[str, set], knockin_loci: set):
    mutation_loci = dict()
    for mut in {"+", "-", "*"}:
        mutation_loci[mut] = candidate_loci[mut] | knockin_loci
    return mutation_loci


def transpose_mutation_loci(mutation_loci: dict[str, set[int]], sequence: str) -> list[set[str]]:
    len_sequence = len(sequence)
    mutation_loci_transposed = [set() for _ in range(len_sequence)]
    for mut, idx_mutation in mutation_loci.items():
        for i, loci in enumerate(mutation_loci_transposed):
            if i in idx_mutation:
                loci.add(mut)
    return mutation_loci_transposed


###########################################################
# main
###########################################################


def cache_indels_count(ARGS, is_control: bool = False, is_insertion: bool = False) -> None:
    tempdir, sample_name, control_name, fasta_alleles = (
        ARGS.tempdir,
        ARGS.sample_name,
        ARGS.control_name,
        ARGS.fasta_alleles,
    )

    for allele, sequence in fasta_alleles.items():
        if is_control and Path(tempdir, control_name, "midsv", f"{allele}_{sample_name}.json").exists():
            prefix = f"{allele}_{sample_name}"
        else:
            prefix = allele

        path_mutation_loci = Path(tempdir, control_name if is_control else sample_name, "mutation_loci")
        if Path(path_mutation_loci, f"{prefix}_count.pickle").exists():
            continue

        path_midsv = Path(tempdir, control_name if is_control else sample_name, "midsv", f"{prefix}.json")
        indels_count, indels_normalized = summarize_indels(path_midsv, sequence)
        io.save_pickle(indels_count, Path(path_mutation_loci, f"{prefix}_count.pickle"))
        io.save_pickle(indels_normalized, Path(path_mutation_loci, f"{prefix}_normalized.pickle"))


def extract_mutation_loci(
    sequence: str,
    path_indels_normalized_sample: Path,
    path_indels_normalized_control: Path,
    path_knockin: Path,
    thresholds: dict[str, float] = {"*": 0.05, "-": 0.05, "+": 0.05},
) -> list[set[str]]:
    indels_normalized_sample = io.load_pickle(path_indels_normalized_sample)
    indels_normalized_control = io.load_pickle(path_indels_normalized_control)

    """Extract candidate mutation loci"""
    indels_normalized_minimize_control = minimize_mutation_counts(indels_normalized_control, indels_normalized_sample)
    anomal_loci = extract_anomal_loci(indels_normalized_sample, indels_normalized_minimize_control, thresholds)

    """Extract error loci in homopolymer regions"""
    errors_in_homopolymer = homopolymer_handler.extract_errors(
        sequence, indels_normalized_sample, indels_normalized_control, anomal_loci
    )
    mutation_loci = discard_errors_in_homopolymer(anomal_loci, errors_in_homopolymer)

    """Merge all mutations and knockin loci"""
    if path_knockin.exists():
        knockin_loci = io.load_pickle(path_knockin)
        mutation_loci = add_knockin_loci(mutation_loci, knockin_loci)

    mutation_loci_merged = merge_index_of_consecutive_indel(mutation_loci)
    mutation_loci_transposed = transpose_mutation_loci(mutation_loci_merged, sequence)
    return mutation_loci_transposed


def cache_mutation_loci(ARGS, is_control: bool = False, is_insertion: bool = False) -> None:
    cache_indels_count(ARGS, is_control, is_insertion)

    if is_control:
        return

    tempdir, sample_name, control_name, fasta_alleles = (
        ARGS.tempdir,
        ARGS.sample_name,
        ARGS.control_name,
        ARGS.fasta_alleles,
    )

    path_mutation_sample = Path(tempdir, sample_name, "mutation_loci")
    path_mutation_control = Path(tempdir, control_name, "mutation_loci")

    for allele, sequence in fasta_alleles.items():
        path_output = Path(path_mutation_sample, f"{allele}.pickle")
        if path_output.exists():
            continue

        file_name = f"{allele}_{sample_name}_normalized.pickle"
        if not Path(path_mutation_control, file_name).exists():
            file_name = f"{allele}_normalized.pickle"
        path_indels_normalized_control = Path(path_mutation_control, file_name)

        path_indels_normalized_sample = Path(path_mutation_sample, f"{allele}_normalized.pickle")
        path_knockin = Path(tempdir, sample_name, "knockin_loci", f"{allele}.pickle")

        mutation_loci = extract_mutation_loci(
            sequence, path_indels_normalized_sample, path_indels_normalized_control, path_knockin
        )

        io.save_pickle(mutation_loci, path_output)
