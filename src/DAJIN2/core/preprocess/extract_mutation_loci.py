from __future__ import annotations

import re
from collections import defaultdict
import json
import numpy as np
from pathlib import Path
from typing import Generator
from scipy import stats
from scipy.spatial import distance
import pickle
from sklearn import linear_model

# from sklearn.neighbors import LocalOutlierFactor
from DAJIN2.core.preprocess.extract_errors_in_homopolymer import extract_errors_in_homopolymer


def read_midsv(filepath) -> Generator[dict[str, str]]:
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
        # Calculate cosine similarity: 1 means exactly same, 0 means completely different.
        # Zero vector does not return correct value, so add 1e-10.
        # example: distance.cosine([0,0,0], [1,2,3]) retuns 0...
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
# Extract dissimilar loci using OneClassSVM
# `extract_different_loci` does not consider the mutation rate in each kmer.
# Thus we got many false positive of kmer with the low percentage of mutation rate
# Consider the mutation rate in the sequence
###########################################################


def _transform_log2(values: np.array) -> np.array:
    values = np.where(values <= 0, 1e-10, values)
    return np.log2(values).reshape(-1, 1)


def _get_divisor(set1, set2) -> int:
    return len(set(set1) & set(set2)) or -1


def _merge_peaks(log2_sample, log2_control, peaks) -> set:
    """Values higher than 75% quantile of the control values and the surrouings are peaks, merge it as a peak"""
    threshold = np.quantile(log2_control, 0.75)
    for i, value in enumerate(log2_sample):
        if i not in peaks and value > threshold:
            for j in range(i - 5, i + 6):
                if j in peaks:
                    peaks.add(i)
                    break
    return peaks


def extract_anomal_loci(indels_sample_normalized, indels_control_normalized) -> dict[str, set]:
    results = dict()
    for mut in ["+", "-", "*"]:
        # preprocess
        values_sample = indels_sample_normalized[mut]
        values_control = indels_control_normalized[mut]
        # values_sample = np.array(indels_sample_normalized[mut])
        # values_control = np.array(indels_control_normalized[mut])
        values_subtract = _transform_log2(values_sample - values_control)
        log2_control = _transform_log2(values_control)
        log2_sample = _transform_log2(values_sample)
        # anomaly detection
        clf = linear_model.SGDOneClassSVM(random_state=0)
        clf.fit(log2_control)
        predicts = clf.predict(log2_sample)
        idx_outliers = np.where(predicts == -1)[0]
        idx_outliers_reverse = np.where(predicts == 1)[0]
        # anomaly detection by quantile
        idx_upper = np.where(values_subtract > np.quantile(log2_control, 0.75))[0]
        # determine which is the correct outliers that is the percentage of outliers is large
        idx1 = len(idx_outliers) / _get_divisor(idx_outliers, idx_upper)
        idx2 = len(idx_outliers_reverse) / _get_divisor(idx_outliers_reverse, idx_upper)
        if idx1 < idx2:
            idx_outliers = idx_outliers_reverse
        results[mut] = _merge_peaks(log2_sample, log2_control, set(idx_outliers) & set(idx_upper))
    return results


###########################################################
# Homolopolymer region
###########################################################


def discard_errors_in_homopolymer(dissimilar_loci, errors_in_homopolymer) -> dict[str, set]:
    mutation_loci = dict()
    for mut in ["+", "-", "*"]:
        error_loci = errors_in_homopolymer[mut]
        mutation_loci[mut] = dissimilar_loci[mut] - error_loci
    return mutation_loci


def transpose_mutation_loci(mutation_loci: set[int], sequence: str) -> list[set]:
    len_sequence = len(sequence)
    mutation_loci_transposed = [set() for _ in range(len_sequence)]
    for mut, idx_mutation in mutation_loci.items():
        for i, loci in enumerate(mutation_loci_transposed):
            if i in idx_mutation:
                loci.add(mut)
    return mutation_loci_transposed


###########################################################
# Biased strand
###########################################################

def discard_errors_on_biased_strand(midsv_sample, mutation_loci) -> dict[str, set]:
    results = dict()
    for mutation, loci in mutation_loci.items():
        mutation_loci_non_biased = set()
        count_plus = defaultdict(int)
        count_total = defaultdict(int)
        midsv_sample = list(midsv_sample)
        for samp in midsv_sample:
            cssplits = samp["CSSPLIT"].split(",")
            for i, cs in enumerate(cssplits):
                if i not in loci:
                    continue
                if not cs.startswith(mutation):
                    continue
                count_total[i] += 1
                if samp["STRAND"] == "+":
                    count_plus[i] += 1
        mutation_loci_non_biased = {i for i, total in count_total.items() if 0.25 < (count_plus[i] / total) < 0.75}
        results[mutation] = mutation_loci_non_biased
    return results

###########################################################
# main
###########################################################


def process_mutation_loci(TEMPDIR: Path, FASTA_ALLELES: dict, CONTROL_NAME: str) -> None:
    for allele, sequence in FASTA_ALLELES.items():
        filepath_control = Path(TEMPDIR, "midsv", f"{CONTROL_NAME}_{allele}.json")
        indels_control = count_indels(read_midsv(filepath_control), sequence)
        coverages_control = call_coverage_on_each_base(read_midsv(filepath_control), sequence)
        indels_control_normalized = normalize_indels(indels_control, coverages_control)
        indels_kmer_control = split_kmer(indels_control_normalized, kmer=11)
        # Save indels_control_normalized and indels_kmer_control as pickle to reuse in consensus calling
        with open(Path(TEMPDIR, "mutation_loci", f"{CONTROL_NAME}_{allele}.pkl"), "wb") as f:
            pickle.dump(indels_control, f)
        with open(Path(TEMPDIR, "mutation_loci", f"{CONTROL_NAME}_{allele}_normalized.pkl"), "wb") as f:
            pickle.dump(indels_control_normalized, f)
        with open(Path(TEMPDIR, "mutation_loci", f"{CONTROL_NAME}_{allele}_kmer.pkl"), "wb") as f:
            pickle.dump(indels_kmer_control, f)


def is_strand_bias(midsv_control) -> bool:
    count_strand = defaultdict(int)
    for m in midsv_control:
        count_strand[m["STRAND"]] += 1
    percentage_plus = count_strand["+"] / (count_strand["+"] + count_strand["-"])
    if 0.25 < percentage_plus < 0.75:
        return False
    else:
        return True


def merge_loci(dissimilar_loci, anomal_loci) -> dict[str, set]:
    mutation_loci = dict()
    for mut in ["+", "-", "*"]:
        mutation_loci[mut] = dissimilar_loci[mut] & anomal_loci[mut]
    return mutation_loci


def extract_mutation_loci(TEMPDIR: Path, FASTA_ALLELES: dict, SAMPLE_NAME: str, CONTROL_NAME: str) -> dict[str, list]:
    MUTATION_LOCI_ALLELES = dict()
    filepath_control = Path(TEMPDIR, "midsv", f"{CONTROL_NAME}_control.json")
    strand_bias = is_strand_bias(read_midsv(filepath_control))
    for allele, sequence in FASTA_ALLELES.items():
        filepath_sample = Path(TEMPDIR, "midsv", f"{SAMPLE_NAME}_{allele}.json")
        indels_sample = count_indels(read_midsv(filepath_sample), sequence)
        coverages_sample = call_coverage_on_each_base(read_midsv(filepath_sample), sequence)
        indels_sample_normalized = normalize_indels(indels_sample, coverages_sample)
        indels_kmer_sample = split_kmer(indels_sample_normalized, kmer=11)
        # Load indels_control_normalized and indels_kmer_control
        with open(Path(TEMPDIR, "mutation_loci", f"{CONTROL_NAME}_{allele}_normalized.pkl"), "rb") as f:
            indels_control_normalized = pickle.load(f)
        with open(Path(TEMPDIR, "mutation_loci", f"{CONTROL_NAME}_{allele}_kmer.pkl"), "rb") as f:
            indels_kmer_control = pickle.load(f)
        # Extract candidate mutation loci
        dissimilar_loci = extract_dissimilar_loci(indels_kmer_sample, indels_kmer_control)
        anomal_loci = extract_anomal_loci(indels_sample_normalized, indels_control_normalized)
        mutation_loci = merge_loci(dissimilar_loci, anomal_loci)
        # Extract error loci in homopolymer regions
        errors_in_homopolymer = extract_errors_in_homopolymer(sequence, indels_sample_normalized, indels_control_normalized, mutation_loci)
        mutation_loci = discard_errors_in_homopolymer(mutation_loci, errors_in_homopolymer)
        if strand_bias is False:
            mutation_loci = discard_errors_on_biased_strand(read_midsv(filepath_sample), mutation_loci)
        mutation_loci_transposed = transpose_mutation_loci(mutation_loci, sequence)
        MUTATION_LOCI_ALLELES[allele] = mutation_loci_transposed
    return MUTATION_LOCI_ALLELES
