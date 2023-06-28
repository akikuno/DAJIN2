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
    extract_anomal_loci,
    merge_loci,
    add_knockin_loci,
    transpose_mutation_loci,
)
from DAJIN2.core.preprocess.extract_errors_in_homopolymer import extract_errors_in_homopolymer


# def _extract_errors_in_homopolymer(indels_sample, indels_control, sequence, candidate_loci) -> dict[str, set]:
#     errors_in_homopolymer = dict()
#     for mut in ["+", "-", "*"]:
#         error_loci: set = extract_errors_in_homopolymer(indels_sample[mut], indels_control[mut], sequence, candidate_loci[mut])
#         errors_in_homopolymer[mut] = error_loci
#     return errors_in_homopolymer


def extract_mutation_loci_by_labels(clust_sample, TEMPDIR, FASTA_ALLELES, CONTROL_NAME, KNOCKIN_LOCI_ALLELES):
    MUTATION_LOCI_LABELS = dict()
    clust_sample.sort(key=lambda x: [x["ALLELE"], x["LABEL"]])
    for (allele, label), group in groupby(clust_sample, key=lambda x: [x["ALLELE"], x["LABEL"]]):
        sequence = FASTA_ALLELES[allele]
        group = list(group)
        coverages_sample = call_coverage_on_each_base(group, sequence)
        indels_sample = count_indels(group, sequence)
        indels_sample_normalized = normalize_indels(indels_sample, coverages_sample)
        indels_kmer_sample = split_kmer(indels_sample_normalized, kmer=11)
        with open(Path(TEMPDIR, "mutation_loci", f"{CONTROL_NAME}_{allele}_normalized.pkl"), "rb") as f:
            indels_control_normalized = pickle.load(f)
        with open(Path(TEMPDIR, "mutation_loci", f"{CONTROL_NAME}_{allele}_kmer.pkl"), "rb") as f:
            indels_kmer_control = pickle.load(f)
        # Calculate dissimilar loci
        dissimilar_loci = extract_dissimilar_loci(indels_kmer_sample, indels_kmer_control)
        anomal_loci = extract_anomal_loci(indels_sample_normalized, indels_control_normalized)
        candidate_loci = merge_loci(dissimilar_loci, anomal_loci)
        # Extract error loci in homopolymer regions
        errors_in_homopolymer = extract_errors_in_homopolymer(sequence, indels_sample_normalized, indels_control_normalized, candidate_loci)
        # errors_in_homopolymer = _extract_errors_in_homopolymer(
        #     indels_sample_normalized, indels_control_normalized, sequence, candidate_loci
        # )
        mutation_loci = discard_errors_in_homopolymer(candidate_loci, errors_in_homopolymer)
        mutation_loci = add_knockin_loci(mutation_loci, KNOCKIN_LOCI_ALLELES[allele])
        mutation_loci_transposed = transpose_mutation_loci(mutation_loci, sequence)
        MUTATION_LOCI_LABELS[label] = mutation_loci_transposed
    return MUTATION_LOCI_LABELS
