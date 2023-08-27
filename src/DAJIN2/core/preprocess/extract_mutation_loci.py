from __future__ import annotations

import re
import pickle
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


def extract_anomal_loci(indels_normalized_sample, indels_normalized_control) -> dict[str, set]:
    results = dict()
    for mut in ["+", "-", "*"]:
        # preprocess
        values_sample = indels_normalized_sample[mut]
        values_control = indels_normalized_control[mut]
        log2_subtract = _transform_log2(values_sample - values_control)
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


def _process_control(TEMPDIR: Path, FASTA_ALLELES: dict, CONTROL_NAME: str) -> None:
    for allele, sequence in FASTA_ALLELES.items():
        path_mutation_loci = Path(TEMPDIR, CONTROL_NAME, "mutation_loci")
        if Path(path_mutation_loci, f"{allele}_count.pickle").exists():
            continue
        filepath_control = Path(TEMPDIR, CONTROL_NAME, "midsv", f"{allele}.json")
        indels_control = count_indels(io.read_jsonl(filepath_control), sequence)
        coverages_control = call_coverage_of_each_base(io.read_jsonl(filepath_control), sequence)
        indels_normalized_control = normalize_indels(indels_control, coverages_control)
        indels_kmer_control = split_kmer(indels_normalized_control, kmer=11)
        # Save indels_normalized_control and indels_kmer_control as pickle to reuse in consensus calling
        with open(Path(path_mutation_loci, f"{allele}_count.pickle"), "wb") as f:
            pickle.dump(indels_control, f)
        with open(Path(path_mutation_loci, f"{allele}_normalized.pickle"), "wb") as f:
            pickle.dump(indels_normalized_control, f)
        with open(Path(path_mutation_loci, f"{allele}_kmer.pickle"), "wb") as f:
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


def extract_mutation_loci(
    TEMPDIR: Path, FASTA_ALLELES: dict, SAMPLE_NAME: str, CONTROL_NAME: str, is_control=False
) -> None:
    if is_control:
        _process_control(TEMPDIR, FASTA_ALLELES, CONTROL_NAME)
        return

    for allele, sequence in FASTA_ALLELES.items():
        path_output = Path(TEMPDIR, SAMPLE_NAME, "mutation_loci", f"{allele}.pickle")
        if path_output.exists():
            continue

        filepath_sample = Path(TEMPDIR, SAMPLE_NAME, "midsv", f"{allele}.json")
        indels_sample = count_indels(io.read_jsonl(filepath_sample), sequence)
        coverages_sample = call_coverage_of_each_base(io.read_jsonl(filepath_sample), sequence)
        indels_normalized_sample = normalize_indels(indels_sample, coverages_sample)
        indels_kmer_sample = split_kmer(indels_normalized_sample, kmer=11)
        # Load indels_normalized_control and indels_kmer_control
        with open(Path(TEMPDIR, CONTROL_NAME, "mutation_loci", f"{allele}_normalized.pickle"), "rb") as f:
            indels_normalized_control = pickle.load(f)
        with open(Path(TEMPDIR, CONTROL_NAME, "mutation_loci", f"{allele}_kmer.pickle"), "rb") as f:
            indels_kmer_control = pickle.load(f)

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
            with open(path_knockin, "rb") as p:
                knockin_loci = pickle.load(p)
            mutation_loci = add_knockin_loci(mutation_loci, knockin_loci)

        mutation_loci = merge_index_of_consecutive_insertions(mutation_loci)
        mutation_loci_transposed = transpose_mutation_loci(mutation_loci, sequence)
        with open(path_output, "wb") as p:
            pickle.dump(mutation_loci_transposed, p)
