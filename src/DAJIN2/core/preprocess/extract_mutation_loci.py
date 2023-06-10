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


def count_indels(midsv_sample, len_sequence: int) -> dict[str, list[int]]:
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


def split_kmer(indels: dict[str, np.array], kmer: int = 10) -> dict[str, np.array]:
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
    Comparing Sample and Control, the 1. 'high similarity',  2. 'similar mean' and
    3. 'similar variance' are considered as sequence errors.
    """
    results = dict()
    for mut in indels_kmer_sample:
        values_sample = indels_kmer_sample[mut]
        values_control = indels_kmer_control[mut]
        # Calculate cosine similarity: 1 means exactly same, 0 means completely different.
        # When calculating cossim, uint32 returns inaccurate results so convert to float64
        cossim = [1 - distance.cosine(x, y) for x, y in zip(values_sample, values_control)]
        # Perform T-test: nan means exactly same, p > 0.05 means similar in average.
        t_pvalues = [stats.ttest_ind(x, y, equal_var=False)[1] for x, y in zip(values_sample, values_control)]
        t_pvalues = [1 if np.isnan(t) else t for t in t_pvalues]
        # Perform F-test: p > 0.05 means similar in variance.
        f_pvalues = [stats.bartlett(x, y)[1] for x, y in zip(values_sample, values_control)]
        # if pvalue == nan or pval > 0.05, samples and controls are similar.
        dissimilar_loci = set()
        for i, (sim, t_pval, f_pval) in enumerate(zip(cossim, t_pvalues, f_pvalues)):
            if not (sim > 0.90 and t_pval > 0.05 and f_pval > 0.05):
                dissimilar_loci.add(i)
        results[mut] = dissimilar_loci
    return results


def discard_errors_in_homopolymer(dissimilar_loci, errors_in_homopolymer) -> dict[str, set]:
    mutation_loci = dict()
    for mut in ["+", "-", "*"]:
        error_loci = errors_in_homopolymer[mut]
        mutation_loci[mut] = dissimilar_loci[mut] - error_loci
    return mutation_loci


def transpose_mutation_loci(mutation_loci: set[int], len_sequence: int) -> list[set]:
    mutation_loci_transposed = [set() for _ in range(len_sequence)]
    for mut, idx_mutation in mutation_loci.items():
        for i, loci in enumerate(mutation_loci_transposed):
            if i in idx_mutation:
                loci.add(mut)
    return mutation_loci_transposed


###########################################################
# main
###########################################################


def extract_mutation_loci(
    TEMPDIR: Path, FASTA_ALLELES: dict, SAMPLE_NAME: str, CONTROL_NAME: str
) -> dict[str, list[set[str]]]:
    MUTATION_LOCI_ALLELES = dict()
    for allele, sequence in FASTA_ALLELES.items():
        filepath_sample = Path(TEMPDIR, "midsv", f"{SAMPLE_NAME}_{allele}.json")
        filepath_control = Path(TEMPDIR, "midsv", f"{CONTROL_NAME}_{allele}.json")
        indels_sample = count_indels(read_midsv(filepath_sample), len(sequence))
        indels_control = count_indels(read_midsv(filepath_control), len(sequence))
        coverages_sample = call_coverage_on_each_base(read_midsv(filepath_sample), sequence)
        coverages_control = call_coverage_on_each_base(read_midsv(filepath_control), sequence)
        indels_sample_normalized = normalize_indels(indels_sample, coverages_sample)
        indels_control_normalized = normalize_indels(indels_control, coverages_control)
        indels_kmer_sample = split_kmer(indels_sample_normalized, kmer=10)
        indels_kmer_control = split_kmer(indels_control_normalized, kmer=10)
        # anomaly_loci = _extract_anomaly_loci(indels_kmer_sample, indels_kmer_control)
        dissimilar_loci = extract_dissimilar_loci(indels_kmer_sample, indels_kmer_control)
        # Extract error loci in homopolymer regions
        errors_in_homopolymer = dict()
        for mut in ["+", "-", "*"]:
            indels_sample_mut = indels_sample_normalized[mut]
            indels_control_mut = indels_control_normalized[mut]
            candidate_loci = dissimilar_loci[mut]
            # candidate_loci = anomaly_loci[mut] & dissimilar_loci[mut]
            errors_in_homopolymer[mut] = extract_errors_in_homopolymer(
                indels_sample_mut, indels_control_mut, sequence, candidate_loci
            )
        mutation_loci = discard_errors_in_homopolymer(dissimilar_loci, errors_in_homopolymer)
        mutation_loci_transposed = transpose_mutation_loci(mutation_loci, len(sequence))
        MUTATION_LOCI_ALLELES[allele] = mutation_loci_transposed
        # Save indels_control_normalized and indels_kmer_control as pickle to reuse in consensus calling
        with open(Path(TEMPDIR, "mutation_loci", f"{CONTROL_NAME}_{allele}_normalized.pkl"), "wb") as f:
            pickle.dump(indels_control_normalized, f)
        with open(Path(TEMPDIR, "mutation_loci", f"{CONTROL_NAME}_{allele}_kmer.pkl"), "wb") as f:
            pickle.dump(indels_kmer_control, f)
    return MUTATION_LOCI_ALLELES
