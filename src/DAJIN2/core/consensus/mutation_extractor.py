from __future__ import annotations

from itertools import groupby
from pathlib import Path

import numpy as np
from sklearn.cluster import MiniBatchKMeans

from DAJIN2.core.consensus.similarity_searcher import cache_selected_control_by_similarity
from DAJIN2.core.preprocess.mutation_extractor import extract_mutation_loci, minimize_mutation_counts, summarize_indels
from DAJIN2.utils import io

"""
Most of the code reuses `preprocess.cache_mutation_loci`.
This code differs from `preprocess.cache_mutation_loci` in that it dynamically extracts control sequences similar to alleles after clustering, and sets thresholds for distinguishing between sequence errors and mutations.
- Extracts reads similar to alleles after clustering to minimize the impact of anomalies in the control sequences.
- Since alleles after clustering essentially assume a single allele, the threshold setting in preprocess is not 10% but is set to a higher value.
"""


def get_thresholds(path_indels_normalized_sample, path_indels_normalized_control) -> dict[str, float]:
    """
    For the consensus threshold, the minimum threshold for "sample - control" was set to be at least 10%.
    """
    indels_normalized_sample = io.load_pickle(path_indels_normalized_sample)
    indels_normalized_control = io.load_pickle(path_indels_normalized_control)
    indels_normalized_minimize_control = minimize_mutation_counts(indels_normalized_control, indels_normalized_sample)
    thresholds = {}
    for mut in {"+", "-", "*"}:
        values_sample = indels_normalized_sample[mut]
        values_control = indels_normalized_minimize_control[mut]
        values_subtract = values_sample - values_control
        kmeans = MiniBatchKMeans(n_clusters=2, random_state=0).fit(values_subtract.reshape(-1, 1))
        threshold = kmeans.cluster_centers_.mean()
        thresholds[mut] = max(threshold, 10)
    return thresholds


def extract_path_control(tempdir: Path, control_name: str, sample_name: str, allele: str) -> Path:
    # Define potential file paths
    paths = [
        Path(tempdir, control_name, "midsv", allele, f"{control_name}.jsonl"),
        Path(tempdir, control_name, "midsv", allele, f"{sample_name}.jsonl"),
    ]
    # Return the first path that exists
    for path in paths:
        if path.exists():
            return path


def extract_path_n_filtered_control(
    tempdir: Path, control_name: str, sample_name: str, path_control: Path, allele: str
) -> Path:
    """
    Filter out reads with N enriched in the control.
    """
    path_output = Path(tempdir, control_name, "consensus", allele, f"{sample_name}_n_filtered.jsonl")
    if path_output.exists():
        return path_output

    path_output.parent.mkdir(parents=True, exist_ok=True)

    midsv_control = io.read_jsonl(path_control)
    n_counts = np.array([sum(1 if cs == "N" else 0 for cs in c["CSSPLIT"].split(",")) for c in midsv_control])

    kmeans = MiniBatchKMeans(n_clusters=2, random_state=0).fit(n_counts.reshape(-1, 1))
    threshold = kmeans.cluster_centers_.mean()
    labels = np.where(n_counts <= threshold, True, False)
    midsv_filtered_control = (m for m, label in zip(io.read_jsonl(path_control), labels) if label)
    io.write_jsonl(midsv_filtered_control, path_output)

    return path_output


def cache_normalized_indels(ARGS, path_consensus_sample: Path, allele: str) -> None:
    sequence = ARGS.fasta_alleles[allele]

    path_control = extract_path_control(ARGS.tempdir, ARGS.control_name, ARGS.sample_name, allele)
    path_control_filtered = extract_path_n_filtered_control(
        ARGS.tempdir, ARGS.control_name, ARGS.sample_name, path_control, allele
    )

    cache_selected_control_by_similarity(ARGS, path_control_filtered, path_consensus_sample, allele)

    path_midsv_similar_control = Path(path_consensus_sample.parent, "control.jsonl")

    _, indels_normalized_sample = summarize_indels(path_consensus_sample, sequence)
    _, indels_normalized_control = summarize_indels(path_midsv_similar_control, sequence)

    io.save_pickle(indels_normalized_sample, Path(path_consensus_sample.parent, "normalized_sample.pickle"))
    io.save_pickle(indels_normalized_control, Path(path_consensus_sample.parent, "normalized_control.pickle"))


def cache_mutation_loci(ARGS, clust_sample: list[dict]) -> None:
    # Separate clusters by label and cache them
    clust_sample.sort(key=lambda x: [x["ALLELE"], x["LABEL"]])
    for (allele, label), group in groupby(clust_sample, key=lambda x: [x["ALLELE"], x["LABEL"]]):
        path_consensus = Path(ARGS.tempdir, ARGS.sample_name, "consensus", allele, str(label))
        path_consensus.mkdir(parents=True, exist_ok=True)

        path_consensus_sample = Path(path_consensus, f"{ARGS.sample_name}_sample.jsonl")
        io.write_jsonl(group, path_consensus_sample)

        cache_normalized_indels(ARGS, path_consensus_sample, allele)

        # Extract and cache mutation loci
        path_indels_normalized_sample = Path(path_consensus, "normalized_sample.pickle")
        path_indels_normalized_control = Path(path_consensus, "normalized_control.pickle")
        sequence = ARGS.fasta_alleles[allele]
        path_knockin = Path(ARGS.tempdir, ARGS.sample_name, "knockin_loci", allele, "knockin.pickle")

        thresholds = get_thresholds(path_indels_normalized_sample, path_indels_normalized_control)

        is_consensus = True

        mutation_loci = extract_mutation_loci(
            sequence,
            path_indels_normalized_sample,
            path_indels_normalized_control,
            path_knockin,
            thresholds,
            is_consensus,
        )

        io.save_pickle(mutation_loci, Path(path_consensus, "mutation_loci.pickle"))
