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
    merge_index_of_consecutive_insertions,
)
from DAJIN2.core.preprocess.extract_errors_in_homopolymer import extract_errors_in_homopolymer


def extract_mutation_loci_by_labels(clust_sample, TEMPDIR, FASTA_ALLELES, CONTROL_NAME, SAMPLE_NAME):
    MUTATION_LOCI_LABELS = dict()
    clust_sample.sort(key=lambda x: [x["ALLELE"], x["LABEL"]])
    for (allele, label), group in groupby(clust_sample, key=lambda x: [x["ALLELE"], x["LABEL"]]):
        sequence = FASTA_ALLELES[allele]
        group = list(group)
        coverages_sample = call_coverage_on_each_base(group, sequence)
        indels_sample = count_indels(group, sequence)
        indels_sample_normalized = normalize_indels(indels_sample, coverages_sample)
        indels_kmer_sample = split_kmer(indels_sample_normalized, kmer=11)
        with open(Path(TEMPDIR, CONTROL_NAME, "mutation_loci", f"{allele}_normalized.pickle"), "rb") as f:
            indels_control_normalized = pickle.load(f)
        with open(Path(TEMPDIR, CONTROL_NAME, "mutation_loci", f"{allele}_kmer.pickle"), "rb") as f:
            indels_kmer_control = pickle.load(f)
        # Calculate dissimilar loci
        dissimilar_loci = extract_dissimilar_loci(indels_kmer_sample, indels_kmer_control)
        anomal_loci = extract_anomal_loci(indels_sample_normalized, indels_control_normalized)
        candidate_loci = merge_loci(dissimilar_loci, anomal_loci)
        # Extract error loci in homopolymer regions
        errors_in_homopolymer = extract_errors_in_homopolymer(
            sequence, indels_sample_normalized, indels_control_normalized, candidate_loci
        )
        mutation_loci = discard_errors_in_homopolymer(candidate_loci, errors_in_homopolymer)
        # Add all mutations into knockin loci
        path_knockin = Path(TEMPDIR, SAMPLE_NAME, "knockin_loci", f"{allele}.pickle")
        if path_knockin.exists():
            with open(path_knockin, "rb") as p:
                knockin_loci = pickle.load(p)
            mutation_loci = add_knockin_loci(mutation_loci, knockin_loci)
        mutation_loci = merge_index_of_consecutive_insertions(mutation_loci)
        mutation_loci_transposed = transpose_mutation_loci(mutation_loci, sequence)
        MUTATION_LOCI_LABELS[label] = mutation_loci_transposed
    return MUTATION_LOCI_LABELS
