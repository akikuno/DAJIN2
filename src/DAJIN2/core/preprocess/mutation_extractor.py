from __future__ import annotations

import re
import numpy as np
from pathlib import Path
from typing import Generator
from collections import defaultdict

from scipy import stats
from scipy.spatial import distance
from sklearn import linear_model

from DAJIN2.utils import io
from DAJIN2.core.preprocess import homopolymer_handler


def call_coverage_of_each_base(midsv_sample: Generator[dict], sequence: str) -> list[int]:
    coverages = [1] * len(sequence)
    for cont in midsv_sample:
        cssplits = cont["CSSPLIT"].split(",")
        for i, cssplit in enumerate(cssplits):
            if cssplit == "N":
                continue
            coverages[i] += 1
    return coverages


def count_indels(midsv_sample, sequence: str) -> dict[str, list[int]]:
    len_sequence = len(sequence)
    count = {"+": [0] * len_sequence, "-": [0] * len_sequence, "*": [0] * len_sequence}
    for samp in midsv_sample:
        for i, cs in enumerate(samp["CSSPLIT"].split(",")):
            if cs.startswith("=") or cs == "N" or re.search(r"a|c|g|t|n", cs):
                continue
            if cs.startswith("+"):
                count["+"][i] += 1
            elif cs.startswith("-"):
                count["-"][i] += 1
            elif cs.startswith("*"):
                count["*"][i] += 1
    return count


def normalize_indels(count: dict[str, list[int]], coverages: list[int]) -> dict[str, np.array]:
    count_normalized = dict()
    coverages = np.array(coverages)
    for mut in count:
        counts = np.array(count[mut])
        count_normalized[mut] = counts / coverages
    return count_normalized


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
# Using Cosine similarity and T test to extract dissimilar Loci
###########################################################


def calculate_cosine_similarities(values_sample: list[float], values_control: list[float]) -> list[float]:
    """
    Calculate cosine similarities between sample and control values.

    Due to the behavior of distance.cosine, when dealing with zero-vectors,
    it doesn't return the expected cosine distance of 1. For example, distance.cosine([0,0,0], [1,2,3]) returns 0.
    To handle this, a small value (1e-10) is added to each element of the vector to prevent them from being zero-vectors.
    This ensures the correct behavior without significantly affecting the cosine similarity calculation.
    """
    return [1 - distance.cosine(x + 1e-10, y + 1e-10) for x, y in zip(values_sample, values_control)]


def perform_t_tests(values_sample: list[float], values_control: list[float]) -> list[float]:
    """
    Perform T-tests between sample and control values.

    If the variance of the samples or control is zero, the p-value of the t-test is NaN.
    In this function, we replace such NaN values with 1, implying that the two samples are similar
    (since a p-value of 1 indicates no statistical difference).
    """
    t_pvalues = [stats.ttest_ind(x, y, equal_var=False)[1] for x, y in zip(values_sample, values_control)]
    return [1 if np.isnan(t) else t for t in t_pvalues]


def find_dissimilar_indices(cossims: list[float], t_pvalues: list[float]) -> set[int]:
    """Identify indices that are dissimilar based on cosine similarities and t-test p-values."""
    return {
        i
        for i, (cossim, t_pval) in enumerate(zip(cossims, t_pvalues))
        if (cossim >= 0.8 and t_pval < 0.05) or cossim < 0.8
    }


def extract_dissimilar_loci(
    indels_kmer_sample: dict[str, list[float]], indels_kmer_control: dict[str, list[float]]
) -> dict[str, set[int]]:
    """
    Compare Sample and Control to identify dissimilar loci.

    Loci that do not closely resemble the reference in both mean and variance, indicating statistically significant differences, are detected as dissimilar loci.
    """
    results = {}
    for mut in {"+", "-", "*"}:
        values_sample = indels_kmer_sample[mut]
        values_control = indels_kmer_control[mut]

        cossims = calculate_cosine_similarities(values_sample, values_control)
        t_pvalues = perform_t_tests(values_sample, values_control)

        results[mut] = find_dissimilar_indices(cossims, t_pvalues)

    return results


###########################################################
# Using OneClassSVM to Extract Anomalous Loci
# The function `extract_dissimilar_loci` overlooks the mutation rate within each kmer.
# As a result, we encounter numerous false positives, especially in kmers with an extremely low mutation rate.
# It's essential to account for the mutation rate across the entire sequence.
###########################################################


def transform_log2(values: np.array) -> np.array:
    """Transform values to log2 scale after handling zeros."""
    values = np.where(values <= 0, 1e-10, values)
    return np.log2(values).reshape(-1, 1)


def merge_surrounding_index(idx_outliers: list[int]) -> set[int]:
    """If an outlier is found in an adjacent 5-mer, the area is also judged as an outlier."""
    idx_merged = set()
    for i, idx_curr in enumerate(idx_outliers):
        if i + 1 == len(idx_outliers):
            break
        idx_next = idx_outliers[i + 1]
        if idx_next - idx_curr <= 5:
            for j in range(idx_curr, idx_next + 1):
                idx_merged.add(j)
        else:
            idx_merged.add(idx_curr)
    return idx_merged


def detect_anomalies(log2_subtract: np.array) -> list[int]:
    """
    Detect anomalies using OneClassSVM and return indices of outliers.

    OneClassSVM classifies data points into two classes: 1 (inliers) and -1 (outliers).
    However, depending on how the "normal" class was learned, either of these classes
    might represent the true anomalies in the context of this problem.

    This function returns the indices of the class with the higher mean of log2_subtract
    values, as this class is considered to be the true anomalies.
    """
    clf = linear_model.SGDOneClassSVM(random_state=0)
    predicts = clf.fit_predict(log2_subtract)
    p1 = [i for i, p in enumerate(predicts) if p == -1]
    p2 = [i for i, p in enumerate(predicts) if p == 1]
    return p1 if np.mean(log2_subtract[p1]) > np.mean(log2_subtract[p2]) else p2


def extract_anomal_loci(indels_normalized_sample, indels_normalized_control) -> dict[str, set[int]]:
    results = {}
    for mut in {"+", "-", "*"}:
        values_sample = indels_normalized_sample[mut]
        values_control = indels_normalized_control[mut]
        log2_subtract = transform_log2(values_sample - values_control)
        idx_outliers = detect_anomalies(log2_subtract)
        results[mut] = merge_surrounding_index(idx_outliers)
    return results


###########################################################
# Homolopolymer region
###########################################################


def discard_errors_in_homopolymer(loci: dict[str, set[int]], errors: dict[str, set[int]]) -> dict[str, set[int]]:
    """Remove detected errors in homopolymer regions from the candidate loci."""
    return {mut: loci[mut] - errors[mut] for mut in {"+", "-", "*"}}


###########################################################
# Merge non-contiguous insertions
###########################################################


def merge_index_of_consecutive_insertions(mutation_loci: dict[str, set[int]]) -> dict[str, set[int]]:
    """Treat as contiguous insertions if there are insertions within five bases of each other"""
    idx_ins_merged = set()
    idx_ins = sorted(mutation_loci["+"])
    for i in range(len(idx_ins) - 1):
        idx_1 = idx_ins[i]
        idx_2 = idx_ins[i + 1]
        if idx_1 + 5 > idx_2:
            for i in range(idx_1, idx_2 + 1):
                idx_ins_merged.add(i)
        else:
            idx_ins_merged.add(idx_1)
    mutation_loci["+"] = idx_ins_merged
    return mutation_loci


###########################################################
# main
###########################################################


def process_data(tempdir: Path, name: str, allele: str, sequence: str) -> tuple:
    """Returns indels, coverages, normalized indels, and kmer indels."""
    path_midsv = Path(tempdir, name, "midsv", f"{allele}.json")
    indels = count_indels(io.read_jsonl(path_midsv), sequence)
    coverages = call_coverage_of_each_base(io.read_jsonl(path_midsv), sequence)
    indels_normalized = normalize_indels(indels, coverages)
    indels_kmer = split_kmer(indels_normalized, kmer=11)

    return indels, indels_normalized, indels_kmer


def merge_loci(dissimilar_loci: dict[str, set], anomal_loci: dict[str, set]) -> dict[str, set]:
    mutation_loci = dict()
    for mut in {"+", "-", "*"}:
        mutation_loci[mut] = dissimilar_loci[mut] & anomal_loci[mut]
    return mutation_loci


def add_knockin_loci(candidate_loci: dict[str, set], knockin_loci: set):
    mutation_loci = dict()
    for mut in {"+", "-", "*"}:
        mutation_loci[mut] = candidate_loci[mut] | knockin_loci
    return mutation_loci


def transpose_mutation_loci(mutation_loci: set[int], sequence: str) -> list[set]:
    len_sequence = len(sequence)
    mutation_loci_transposed = [set() for _ in range(len_sequence)]
    for mut, idx_mutation in mutation_loci.items():
        for i, loci in enumerate(mutation_loci_transposed):
            if i in idx_mutation:
                loci.add(mut)
    return mutation_loci_transposed


def extract_mutation_loci(
    TEMPDIR: Path, FASTA_ALLELES: dict, SAMPLE_NAME: str, CONTROL_NAME: str, is_control=False
) -> None:
    path_mutation_cont = Path(TEMPDIR, CONTROL_NAME, "mutation_loci")
    for allele, sequence in FASTA_ALLELES.items():
        if is_control:
            if Path(path_mutation_cont, f"{allele}_count.pickle").exists():
                continue
            indels_control, indels_normalized_control, indels_kmer_control = process_data(
                TEMPDIR, CONTROL_NAME, allele, sequence
            )

            # Save control data for later use
            io.save_pickle(indels_control, Path(path_mutation_cont, f"{allele}_count.pickle"))
            io.save_pickle(indels_normalized_control, Path(path_mutation_cont, f"{allele}_normalized.pickle"))
            io.save_pickle(indels_kmer_control, Path(path_mutation_cont, f"{allele}_kmer.pickle"))
            continue

        path_output = Path(TEMPDIR, SAMPLE_NAME, "mutation_loci", f"{allele}.pickle")
        if path_output.exists():
            continue

        _, indels_normalized_sample, indels_kmer_sample = process_data(TEMPDIR, SAMPLE_NAME, allele, sequence)

        # Load indels_normalized_control and indels_kmer_control
        indels_normalized_control = io.load_pickle(Path(path_mutation_cont, f"{allele}_normalized.pickle"))
        indels_kmer_control = io.load_pickle(Path(path_mutation_cont, f"{allele}_kmer.pickle"))

        # Extract candidate mutation loci
        dissimilar_loci = extract_dissimilar_loci(indels_kmer_sample, indels_kmer_control)
        anomal_loci = extract_anomal_loci(indels_normalized_sample, indels_normalized_control)
        candidate_loci = merge_loci(dissimilar_loci, anomal_loci)

        # Extract error loci in homopolymer regions
        errors_in_homopolymer = homopolymer_handler.extract_errors(
            sequence, indels_normalized_sample, indels_normalized_control, candidate_loci
        )
        mutation_loci = discard_errors_in_homopolymer(candidate_loci, errors_in_homopolymer)

        # Merge all mutations and knockin loci
        path_knockin = Path(TEMPDIR, SAMPLE_NAME, "knockin_loci", f"{allele}.pickle")
        if path_knockin.exists():
            knockin_loci = io.load_pickle(path_knockin)
            mutation_loci = add_knockin_loci(mutation_loci, knockin_loci)

        mutation_loci = merge_index_of_consecutive_insertions(mutation_loci)
        mutation_loci_transposed = transpose_mutation_loci(mutation_loci, sequence)
        io.save_pickle(mutation_loci_transposed, path_output)
