from __future__ import annotations

from pathlib import Path
from itertools import groupby
import pickle

from DAJIN2.core.preprocess.extract_mutation_loci import (
    discard_errors_in_homopolymer,
    call_coverage_on_each_base,
    count_indels,
    normalize_indels,
    split_kmer,
    extract_dissimilar_loci,
    transpose_mutation_loci,
)
from DAJIN2.core.preprocess.extract_errors_in_homopolymer import extract_errors_in_homopolymer


def _extract_errors_in_homopolymer(indels_sample, indels_control, sequence, dissimilar_loci) -> dict[str, set]:
    errors_in_homopolymer = dict()
    for mut in ["+", "-", "*"]:
        indels_sample_mut = indels_sample[mut]
        indels_control_mut = indels_control[mut]
        candidate_loci = dissimilar_loci[mut]
        error_loci: set = extract_errors_in_homopolymer(indels_sample_mut, indels_control_mut, sequence, candidate_loci)
        errors_in_homopolymer[mut] = error_loci
    return errors_in_homopolymer


def extract_mutation_loci_by_labels(clust_sample, TEMPDIR, FASTA_ALLELES, CONTROL_NAME):
    MUTATION_LOCI_LABELS = dict()
    clust_sample.sort(key=lambda x: [x["ALLELE"], x["LABEL"]])
    for (allele, label), group in groupby(clust_sample, key=lambda x: [x["ALLELE"], x["LABEL"]]):
        sequence = FASTA_ALLELES[allele]
        group = list(group)
        coverages_sample = call_coverage_on_each_base(group, sequence)
        indels_sample = count_indels(group, sequence)
        indels_sample_normalized = normalize_indels(indels_sample, coverages_sample)
        indels_kmer_sample = split_kmer(indels_sample_normalized, kmer=10)
        with open(Path(TEMPDIR, "mutation_loci", f"{CONTROL_NAME}_{allele}_normalized.pkl"), "rb") as f:
            indels_control_normalized = pickle.load(f)
        with open(Path(TEMPDIR, "mutation_loci", f"{CONTROL_NAME}_{allele}_kmer.pkl"), "rb") as f:
            indels_kmer_control = pickle.load(f)
        dissimilar_loci: dict = extract_dissimilar_loci(indels_kmer_sample, indels_kmer_control)
        # Extract error loci in homopolymer regions
        errors_in_homopolymer = _extract_errors_in_homopolymer(
            indels_sample_normalized, indels_control_normalized, sequence, dissimilar_loci
        )
        mutation_loci = discard_errors_in_homopolymer(dissimilar_loci, errors_in_homopolymer)
        mutation_loci_transposed = transpose_mutation_loci(mutation_loci, sequence)
        MUTATION_LOCI_LABELS[label] = mutation_loci_transposed
    return MUTATION_LOCI_LABELS
