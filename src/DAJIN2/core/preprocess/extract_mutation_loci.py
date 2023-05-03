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


def _count_indels(midsv_sample, len_sequence: int) -> dict[str, list[int]]:
    count = {"+": [0] * len_sequence,
            "-": [0] * len_sequence,
            "*": [0] * len_sequence}
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


def _split_kmer(indels: dict[str, list[int]], kmer: int = 10) -> dict[str, list[list[int]]]:
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
                results[mut].append(value[start : end])
            else:
                results[mut].append([0]*kmer)
    return results


def _extract_anomaly_loci(indels_kmer_sample: dict, indels_kmer_control: dict, coverage_sample: int, coverage_control: int) -> dict[str, set[int]]:
    anomaly_loci = dict()
    clf = LocalOutlierFactor(novelty=True, n_neighbors=5)
    for key in indels_kmer_sample.keys():
        loci = set()
        values_control = np.array(indels_kmer_control[key]) / coverage_control
        values_sample = np.array(indels_kmer_sample[key]) / coverage_sample
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
        anomaly_loci.update({key: loci})
    return anomaly_loci


def _extract_dissimilar_loci(indels_kmer_sample: dict[str, list[list[int]]], indels_kmer_control: dict[str, list[list[int]]]) -> dict[str, set]:
    results = dict()
    for mut in indels_kmer_sample:
        cossim = [distance.cosine(x, y) for x, y in zip(indels_kmer_sample[mut], indels_kmer_control[mut])]
        pvalues = [stats.ttest_ind(x, y, equal_var=False)[1] for x, y in zip(indels_kmer_sample[mut], indels_kmer_control[mut])]
        # if pvalue == nan, samples and controls are exactly same.
        cossim_pval_false = [cossim if pvalue > 0.05 or np.isnan(pvalue) else 1 for cossim, pvalue in zip(cossim, pvalues)]
        dissimilar_loci = {i for i, x in enumerate(cossim_pval_false) if x > 0.05}
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


def read_midsv(filepath) -> Generator[dict[str, str]]:
    with open(filepath, 'r') as f:
        for line in f:
            yield json.loads(line)


def count_newlines(filepath):
    def _make_gen(reader):
        while True:
            b = reader(2 ** 16)
            if not b: break
            yield b
    with open(filepath, "rb") as f:
        count = sum(buf.count(b"\n") for buf in _make_gen(f.raw.read))
    return count


def extract_mutation_loci(TEMPDIR: Path, FASTA_ALLELES: dict, SAMPLE_NAME: str, CONTROL_NAME:str) -> dict[str, list[set[str]]]:
    MUTATION_LOCI_ALLELES = dict()
    for allele, sequence in FASTA_ALLELES.items():
        filepath_sample = Path(TEMPDIR, "midsv", f"{SAMPLE_NAME}_{allele}.json")
        filepath_control = Path(TEMPDIR, "midsv", f"{CONTROL_NAME}_{allele}.json")
        covarage_sample = count_newlines(filepath_sample)
        covarage_control = count_newlines(filepath_control)
        indels_sample = _count_indels(read_midsv(filepath_sample), len(sequence))
        indels_control = _count_indels(read_midsv(filepath_control), len(sequence))
        indels_kmer_sample = _split_kmer(indels_sample, kmer = 10)
        indels_kmer_control = _split_kmer(indels_control, kmer = 10)
        anomaly_loci = _extract_anomaly_loci(indels_kmer_sample, indels_kmer_control, covarage_sample, covarage_control)
        dissimilar_loci = _extract_dissimilar_loci(indels_kmer_sample, indels_kmer_control)
        mutation_loci = dict()
        for mut in anomaly_loci:
            mutation_loci.update({mut: anomaly_loci[mut] & dissimilar_loci[mut]})
        mutation_loci_transposed = _transpose_mutation_loci(mutation_loci, len(sequence))
        MUTATION_LOCI_ALLELES.update({allele: mutation_loci_transposed})
    return MUTATION_LOCI_ALLELES
