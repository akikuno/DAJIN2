from __future__ import annotations

import re
import json
import pickle
import numpy as np
from pathlib import Path
from typing import Generator
from collections import defaultdict
from scipy import stats
from scipy.spatial import distance
from sklearn import linear_model

from DAJIN2.core.preprocess.extract_errors_in_homopolymer import extract_errors_in_homopolymer


def read_midsv(filepath: str | Path) -> Generator[dict[str, str]]:
    with open(filepath, "r") as f:
        for line in f:
            yield json.loads(line)


def call_coverage_on_each_base(midsv_sample: Generator[dict], sequence: str) -> list[int]:
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
                # count["+"][i] += len(cs.split("|"))
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


def extract_dissimilar_loci(indels_kmer_sample: dict, indels_kmer_control: dict) -> dict[str, set]:
    """
    Comparing Sample and Control, the 'similar mean' and
    'similar variance' are considered as sequence errors.
    """
    results = dict()
    for mut in ["+", "-", "*"]:
        values_sample = indels_kmer_sample[mut]
        values_control = indels_kmer_control[mut]
        """Calculate cosine similarity: 1 means exactly same, 0 means completely different.
        - Zero vector does not return correct value, so add 1e-10.
            - example: distance.cosine([0,0,0], [1,2,3]) == 0
        """
        cossims = [1 - distance.cosine(x + 1e-10, y + 1e-10) for x, y in zip(values_sample, values_control)]
        # Perform T-test: nan means exactly same, p > 0.05 means similar in average.
        t_pvalues = [stats.ttest_ind(x, y, equal_var=False)[1] for x, y in zip(values_sample, values_control)]
        t_pvalues = [1 if np.isnan(t) else t for t in t_pvalues]
        # if pvalue == nan or pval > 0.05, samples and controls are similar.
        dissimilar_loci = set()
        for i, (cossim, t_pval) in enumerate(zip(cossims, t_pvalues)):
            if (cossim >= 0.8 and t_pval < 0.05) or cossim < 0.8:
                dissimilar_loci.add(i)
        results[mut] = dissimilar_loci
    return results


###########################################################
# Extract anomaly loci using OneClassSVM
# `extract_different_loci` does not consider the mutation rate in each kmer.
# Thus we got many false positives of kmer with the very low percentage of mutation rate
# Consider the mutation rate in the whole sequence
###########################################################


def _transform_log2(values: np.array) -> np.array:
    values = np.where(values <= 0, 1e-10, values)
    return np.log2(values).reshape(-1, 1)


# def _ratio_of_outliers(idx_outliers, idx_upper) -> int:
#     """
#     Determine which is the correct outliers (1 or -1) by the ratio of outliers
#     """
#     if len(idx_outliers) == 0:
#         return -1
#     return len(set(idx_outliers) & set(idx_upper)) / len(idx_outliers)


# def _get_divisor(set1, set2) -> int:
#     return len(set(set1) & set(set2)) or -1


def _merge_surrounding_index(idx_outliers: list) -> set:
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


# def _merge_peaks(log2_sample, log2_control, peaks) -> set:
#     """Values higher than 75% quantile of the control values and
#     the surrouings are peaks, merge it as a peak"""
#     threshold = np.quantile(log2_control, 0.75)
#     for i, value in enumerate(log2_sample):
#         if i not in peaks and value > threshold:
#             for j in range(i - 5, i + 6):
#                 if j in peaks:
#                     peaks.add(i)
#                     break
#     return peaks


def extract_anomal_loci(indels_sample_normalized, indels_control_normalized) -> dict[str, set]:
    results = dict()
    for mut in ["+", "-", "*"]:
        # preprocess
        values_sample = indels_sample_normalized[mut]
        values_control = indels_control_normalized[mut]
        log2_subtract = _transform_log2(values_sample - values_control)
        # log2_control = _transform_log2(values_control)
        # log2_sample = _transform_log2(values_sample)
        # anomaly detection
        clf = linear_model.SGDOneClassSVM(random_state=0)
        predicts = clf.fit_predict(log2_subtract)
        # clf.fit(log2_control)
        # predicts = clf.predict(log2_sample)
        p1 = [i for i, p in enumerate(predicts) if p == -1]
        p2 = [i for i, p in enumerate(predicts) if p == 1]
        if np.mean(log2_subtract[p1]) > np.mean(log2_subtract[p2]):
            idx_outliers = p1
        else:
            idx_outliers = p2
        results[mut] = _merge_surrounding_index(idx_outliers)
        # # anomaly detection by quantile
        # threshold = np.median(np.concatenate([log2_control[:100], log2_control[-100:]]))
        # idx_upper = [i for i, s in enumerate(log2_sample) if s > threshold]
        # # determine which is the correct outliers that is the percentage of outliers is large
        # if _ratio_of_outliers(idx_outliers, idx_upper) < _ratio_of_outliers(idx_outliers_reverse, idx_upper):
        #     idx_outliers = idx_outliers_reverse
        # results[mut] = _merge_peaks(log2_sample, log2_control, set(idx_outliers) & set(idx_upper))
    return results


###########################################################
# Homolopolymer region
###########################################################


def discard_errors_in_homopolymer(candidate_loci: dict[str, set], errors: dict[str, set]) -> dict[str, set]:
    mutation_loci = dict()
    for mut in ["+", "-", "*"]:
        mutation_loci[mut] = candidate_loci[mut] - errors[mut]
    return mutation_loci


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


def process_mutation_loci(TEMPDIR: Path, FASTA_ALLELES: dict, CONTROL_NAME: str) -> None:
    for allele, sequence in FASTA_ALLELES.items():
        if Path(TEMPDIR, "mutation_loci", f"{CONTROL_NAME}_{allele}.pickle").exists():
            continue
        filepath_control = Path(TEMPDIR, "midsv", f"{CONTROL_NAME}_{allele}.json")
        indels_control = count_indels(read_midsv(filepath_control), sequence)
        coverages_control = call_coverage_on_each_base(read_midsv(filepath_control), sequence)
        indels_control_normalized = normalize_indels(indels_control, coverages_control)
        indels_kmer_control = split_kmer(indels_control_normalized, kmer=11)
        # Save indels_control_normalized and indels_kmer_control as pickle to reuse in consensus calling
        with open(Path(TEMPDIR, "mutation_loci", f"{CONTROL_NAME}_{allele}_count.pickle"), "wb") as f:
            pickle.dump(indels_control, f)
        with open(Path(TEMPDIR, "mutation_loci", f"{CONTROL_NAME}_{allele}_normalized.pickle"), "wb") as f:
            pickle.dump(indels_control_normalized, f)
        with open(Path(TEMPDIR, "mutation_loci", f"{CONTROL_NAME}_{allele}_kmer.pickle"), "wb") as f:
            pickle.dump(indels_kmer_control, f)


def merge_loci(dissimilar_loci: dict[str, set], anomal_loci: dict[str, set]) -> dict[str, set]:
    mutation_loci = dict()
    for mut in ["+", "-", "*"]:
        mutation_loci[mut] = dissimilar_loci[mut] & anomal_loci[mut]
    return mutation_loci


def add_knockin_loci(candidate_loci: dict[str, set], knockin_loci: set):
    mutation_loci = dict()
    for mut in ["+", "-", "*"]:
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


def extract_mutation_loci(TEMPDIR: Path, FASTA_ALLELES: dict, SAMPLE_NAME: str, CONTROL_NAME: str) -> None:
    for allele, sequence in FASTA_ALLELES.items():
        if Path(TEMPDIR, "mutation_loci", f"{SAMPLE_NAME}_{allele}.pickle").exists():
            continue
        filepath_sample = Path(TEMPDIR, "midsv", f"{SAMPLE_NAME}_{allele}.json")
        indels_sample = count_indels(read_midsv(filepath_sample), sequence)
        coverages_sample = call_coverage_on_each_base(read_midsv(filepath_sample), sequence)
        indels_sample_normalized = normalize_indels(indels_sample, coverages_sample)
        indels_kmer_sample = split_kmer(indels_sample_normalized, kmer=11)
        # Load indels_control_normalized and indels_kmer_control
        with open(Path(TEMPDIR, "mutation_loci", f"{CONTROL_NAME}_{allele}_normalized.pickle"), "rb") as f:
            indels_control_normalized = pickle.load(f)
        with open(Path(TEMPDIR, "mutation_loci", f"{CONTROL_NAME}_{allele}_kmer.pickle"), "rb") as f:
            indels_kmer_control = pickle.load(f)
        # Extract candidate mutation loci
        dissimilar_loci = extract_dissimilar_loci(indels_kmer_sample, indels_kmer_control)
        anomal_loci = extract_anomal_loci(indels_sample_normalized, indels_control_normalized)
        candidate_loci = merge_loci(dissimilar_loci, anomal_loci)
        # Extract error loci in homopolymer regions
        errors_in_homopolymer = extract_errors_in_homopolymer(
            sequence, indels_sample_normalized, indels_control_normalized, candidate_loci
        )
        mutation_loci = discard_errors_in_homopolymer(candidate_loci, errors_in_homopolymer)
        # Add all mutations into knockin loci
        if Path(TEMPDIR, "knockin_loci", f"{SAMPLE_NAME}_{allele}.pickle").exists():
            with open(Path(TEMPDIR, "knockin_loci", f"{SAMPLE_NAME}_{allele}.pickle"), "rb") as p:
                knockin_loci = pickle.load(p)
            mutation_loci = add_knockin_loci(mutation_loci, knockin_loci)
        mutation_loci = merge_index_of_consecutive_insertions(mutation_loci)
        mutation_loci_transposed = transpose_mutation_loci(mutation_loci, sequence)
        with open(Path(TEMPDIR, "mutation_loci", f"{SAMPLE_NAME}_{allele}.pickle"), "wb") as p:
            pickle.dump(mutation_loci_transposed, p)
