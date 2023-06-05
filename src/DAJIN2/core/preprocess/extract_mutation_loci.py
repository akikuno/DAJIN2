from __future__ import annotations

import re
from collections import defaultdict
import json
import numpy as np
from pathlib import Path
from typing import Generator
from scipy import stats
from scipy.spatial import distance
from sklearn.neighbors import LocalOutlierFactor
from DAJIN2.core.preprocess.extract_errors_in_homopolymer import extract_errors_in_homopolymer


def read_midsv(filepath) -> Generator[dict[str, str]]:
    with open(filepath, "r") as f:
        for line in f:
            yield json.loads(line)


def _make_gen(reader):
    while True:
        b = reader(2**16)
        if not b:
            break
        yield b


def count_newlines(filepath):
    with open(filepath, "rb") as f:
        count = sum(buf.count(b"\n") for buf in _make_gen(f.raw.read))
    return count


def _count_indels(midsv_sample, len_sequence: int) -> dict[str, list[int]]:
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


def _normalize_indels(count: dict[str, list[int]], coverage: int) -> dict[str, np.array]:
    count_normalized = dict()
    for mut in count:
        counts = np.array(count[mut])
        counts_norm = (counts / coverage) * 10**6
        count_normalized[mut] = counts_norm.astype("uint32")
    return count_normalized


def _split_kmer(indels: dict[str, np.array], kmer: int = 10) -> dict[str, np.array]:
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
                results[mut].append(np.array([0] * kmer, dtype="uint32"))
    return results


def _extract_anomaly_loci(indels_kmer_sample: dict, indels_kmer_control: dict) -> dict[str, set[int]]:
    anomaly_loci = dict()
    clf = LocalOutlierFactor(novelty=True, n_neighbors=5)
    for mut in indels_kmer_sample.keys():
        loci = set()
        values_control = indels_kmer_control[mut]
        values_sample = indels_kmer_sample[mut]
        index = -1
        for i, (value_control, value_sample) in enumerate(zip(values_control, values_sample)):
            if i == index:
                continue
            clf.fit(value_control.reshape(-1, 1))
            pred = clf.predict(value_sample.reshape(-1, 1))
            if pred[5] == -1:
                loci.add(i)
            # If the next base is not -1, do not validate the next base because the next base is not an outlier.
            if pred[6] == 1:
                index = i + 1
        anomaly_loci.update({mut: loci})
    return anomaly_loci


def _extract_dissimilar_loci(indels_kmer_sample: dict, indels_kmer_control: dict) -> dict[str, set]:
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
        cossim = [
            1 - distance.cosine(x.astype("float64"), y.astype("float64")) for x, y in zip(values_sample, values_control)
        ]
        # Perform T-test: nan means exactly same, p > 0.05 means similar in average.
        t_pvalues = [stats.ttest_ind(x, y, equal_var=False)[1] for x, y in zip(values_sample, values_control)]
        t_pvalues = [1 if np.isnan(t) else t for t in t_pvalues]
        # Perform F-test: p > 0.05 means similar in variance.
        f_pvalues = [stats.bartlett(x, y)[1] for x, y in zip(values_sample, values_control)]
        # if pvalue == nan or pval > 0.05, samples and controls are similar.
        dissimilar_loci = set()
        for i, (sim, t_pval, f_pval) in enumerate(zip(cossim, t_pvalues, f_pvalues)):
            flag_seqerror = False
            if sim > 0.90 and t_pval > 0.05 and f_pval > 0.05:
                flag_seqerror = True
            if flag_seqerror is False:
                dissimilar_loci.add(i)
        results.update({mut: dissimilar_loci})
    return results


def _transpose_mutation_loci(mutation_loci, len_sequence):
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
        coverage_sample = count_newlines(filepath_sample)
        coverage_control = count_newlines(filepath_control)
        indels_sample = _count_indels(read_midsv(filepath_sample), len(sequence))
        indels_control = _count_indels(read_midsv(filepath_control), len(sequence))
        indels_sample_normalized = _normalize_indels(indels_sample, coverage_sample)
        indels_control_normalized = _normalize_indels(indels_control, coverage_control)
        indels_kmer_sample = _split_kmer(indels_sample_normalized, kmer=10)
        indels_kmer_control = _split_kmer(indels_control_normalized, kmer=10)
        anomaly_loci = _extract_anomaly_loci(indels_kmer_sample, indels_kmer_control)
        dissimilar_loci = _extract_dissimilar_loci(indels_kmer_sample, indels_kmer_control)
        # Extract error loci in homopolymer regions
        error_loci_homopolymer = dict()
        for mut in ["+", "-", "*"]:
            candidate_loci = anomaly_loci[mut] & dissimilar_loci[mut]
            indels_sample_mut = indels_sample[mut]
            indels_control_mut = indels_control[mut]
            error_loci = extract_errors_in_homopolymer(indels_sample_mut, indels_control_mut, sequence, candidate_loci)
            error_loci_homopolymer.update({mut: error_loci})
        mutation_loci = dict()
        for mut in ["+", "-", "*"]:
            candidate_loci = anomaly_loci[mut] & dissimilar_loci[mut]
            error_loci = error_loci_homopolymer[mut]
            # Discard error loci in homopolymer regions
            mutation_loci.update({mut: candidate_loci - error_loci})
        mutation_loci_transposed = _transpose_mutation_loci(mutation_loci, len(sequence))
        MUTATION_LOCI_ALLELES.update({allele: mutation_loci_transposed})
    return MUTATION_LOCI_ALLELES
