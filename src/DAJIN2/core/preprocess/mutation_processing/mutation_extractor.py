from __future__ import annotations

from pathlib import Path

from DAJIN2.core.preprocess.error_correction.homopolymer_handler import extract_sequence_errors_in_homopolymer_loci
from DAJIN2.core.preprocess.error_correction.sequence_error_handler import (
    extract_sequence_errors_using_insample_control,
)
from DAJIN2.core.preprocess.error_correction.strand_bias_handler import extract_sequence_errors_in_strand_biased_loci
from DAJIN2.core.preprocess.mutation_processing.anomaly_detector import extract_anomal_loci
from DAJIN2.core.preprocess.mutation_processing.indel_counter import minimize_mutation_counts, summarize_indels
from DAJIN2.core.preprocess.mutation_processing.indel_merger import (
    add_knockin_loci,
    discard_errors,
    merge_index_of_consecutive_indel,
    transpose_mutation_loci,
)
from DAJIN2.utils import config, io

"""
To suppress the following warnings from `scipy.wilcoxon`:
UserWarning: Exact p-value calculation does not work if there are zeros.
"""
config.set_warnings_ignore()


def extract_mutation_loci(
    path_midsv_sample: Path,
    sequence: str,
    path_indels_normalized_sample: Path,
    path_indels_normalized_control: Path,
    path_knockin: Path,
    thresholds: dict[str, float] = None,
    is_consensus: bool = False,
) -> list[set[str]]:
    """
    Extract mutation loci from sample and control data.

    This is the core mutation extraction function for preprocess module.
    """
    if thresholds is None:
        thresholds = {"*": 0.1, "-": 0.1, "+": 0.1}
    indels_normalized_sample: dict[str, list[float]] = io.load_pickle(path_indels_normalized_sample)

    # Extract candidate mutation loci
    indels_normalized_control: dict[str, list[float]] = minimize_mutation_counts(
        io.load_pickle(path_indels_normalized_control), indels_normalized_sample
    )
    anomal_loci: dict[str, set[int]] = extract_anomal_loci(
        indels_normalized_sample, indels_normalized_control, thresholds, is_consensus
    )

    # Merge all mutations and knockin loci
    if path_knockin.exists():
        knockin_loci = io.load_pickle(path_knockin)
        anomal_loci = add_knockin_loci(anomal_loci, knockin_loci)

    # Extract error loci in homopolymer regions
    errors_in_homopolymer: dict[str, set[int]] = extract_sequence_errors_in_homopolymer_loci(
        sequence, indels_normalized_sample, indels_normalized_control, anomal_loci
    )
    # Extract strand biased loci
    bias_in_strandness: dict[str, set[int]] = extract_sequence_errors_in_strand_biased_loci(
        path_midsv_sample, transpose_mutation_loci(anomal_loci, sequence)
    )

    anomal_loci = discard_errors(anomal_loci, errors_in_homopolymer)
    anomal_loci = discard_errors(anomal_loci, bias_in_strandness)

    if is_consensus:
        # Extract errors using in-sample control
        # In consensus analysis, errors are removed more strictly.
        # Reason not applied to non-consensus analysis:
        # Applying overly strict filtering from the beginning may cause minor variants
        # at around 1% frequency to be mistakenly removed as errors.
        # In contrast, in consensus analysis, if a 1% allele is properly separated,
        # it will be merged as a 100% allele at that stage,
        # and thus should not be treated as an error.
        errors_using_insample_control: dict[str, set[int]] = extract_sequence_errors_using_insample_control(
            indels_normalized_sample, sigma_threshold=2.0
        )
        anomal_loci = discard_errors(anomal_loci, errors_using_insample_control)

    anomal_loci_merged = merge_index_of_consecutive_indel(anomal_loci)
    mutation_loci = transpose_mutation_loci(anomal_loci_merged, sequence)
    return mutation_loci


def cache_indels_count(ARGS, is_control: bool = False) -> None:
    """Cache indel counts for samples or controls."""
    dirname = ARGS.control_name if is_control else ARGS.sample_name
    for allele, sequence in ARGS.fasta_alleles.items():
        path_mutation_loci = Path(ARGS.tempdir, dirname, "mutation_loci", allele)
        path_mutation_loci.mkdir(parents=True, exist_ok=True)

        if not is_control:
            prefix = ARGS.sample_name
        else:
            path_insertion = Path(ARGS.tempdir, ARGS.control_name, "midsv", allele, f"{ARGS.sample_name}.jsonl")
            if path_insertion.exists():
                prefix = ARGS.sample_name
            else:
                prefix = ARGS.control_name

        if Path(path_mutation_loci, f"{prefix}_count.pickle").exists():
            continue

        path_midsv = Path(ARGS.tempdir, dirname, "midsv", allele, f"{prefix}.jsonl")
        indels_count, indels_normalized = summarize_indels(path_midsv, sequence)
        io.save_pickle(indels_count, Path(path_mutation_loci, f"{prefix}_count.pickle"))
        io.save_pickle(indels_normalized, Path(path_mutation_loci, f"{prefix}_normalized.pickle"))


def cache_mutation_loci(ARGS, is_control: bool = False) -> None:
    """Cache mutation loci for preprocess stage."""
    cache_indels_count(ARGS, is_control)

    if is_control:
        return None

    for allele, sequence in ARGS.fasta_alleles.items():
        path_midsv_sample = Path(ARGS.tempdir, ARGS.sample_name, "midsv", allele, f"{ARGS.sample_name}.jsonl")
        path_mutation_sample = Path(ARGS.tempdir, ARGS.sample_name, "mutation_loci", allele)
        path_mutation_control = Path(ARGS.tempdir, ARGS.control_name, "mutation_loci", allele)

        path_output_mutation_loci = Path(path_mutation_sample, "mutation_loci.pickle")
        if path_output_mutation_loci.exists():
            continue

        file_name = f"{ARGS.sample_name}_normalized.pickle"
        if not Path(path_mutation_control, file_name).exists():
            file_name = f"{ARGS.control_name}_normalized.pickle"

        path_indels_normalized_control = Path(path_mutation_control, file_name)
        path_indels_normalized_sample = Path(path_mutation_sample, f"{ARGS.sample_name}_normalized.pickle")
        path_knockin = Path(ARGS.tempdir, ARGS.sample_name, "knockin_loci", allele, "knockin.pickle")

        mutation_loci: list[set[str]] = extract_mutation_loci(
            path_midsv_sample, sequence, path_indels_normalized_sample, path_indels_normalized_control, path_knockin
        )

        io.save_pickle(mutation_loci, path_output_mutation_loci)
