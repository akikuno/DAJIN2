from __future__ import annotations

from pathlib import Path
from itertools import groupby

import numpy as np
from sklearn.cluster import MiniBatchKMeans

from DAJIN2.utils import io
from DAJIN2.core.preprocess.mutation_extractor import summarize_indels, extract_mutation_loci, minimize_mutation_counts
from DAJIN2.core.consensus.similarity_searcher import cache_selected_control_by_similarity

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
    thresholds = dict()
    for mut in {"+", "-", "*"}:
        values_sample = indels_normalized_sample[mut]
        values_control = indels_normalized_minimize_control[mut]
        values_subtract = values_sample - values_control
        kmeans = MiniBatchKMeans(n_clusters=2, random_state=0).fit(values_subtract.reshape(-1, 1))
        threshold = kmeans.cluster_centers_.mean()
        thresholds[mut] = max(threshold, 10)
    return thresholds


def extract_path_control(ARGS, allele: str) -> Path:
    # Define potential file paths
    paths = [
        Path(ARGS.tempdir, ARGS.control_name, "midsv", f"{allele}.json"),
        Path(ARGS.tempdir, ARGS.control_name, "midsv", f"{allele}_{ARGS.sample_name}.json"),
    ]
    # Return the first path that exists
    for path in paths:
        if path.exists():
            return path


def extract_path_n_filtered_control(ARGS, path_midsv_control: Path) -> Path:
    """
    Filter out reads with N enriched in the control.
    """
    path_output = Path(ARGS.tempdir, ARGS.control_name, "consensus", f"{path_midsv_control.stem}_n_filtered.jsonl")
    if path_output.exists():
        return path_output

    midsv_control = io.read_jsonl(path_midsv_control)
    n_counts = np.array([sum(1 if cs == "N" else 0 for cs in c["CSSPLIT"].split(",")) for c in midsv_control])

    kmeans = MiniBatchKMeans(n_clusters=2, random_state=0).fit(n_counts.reshape(-1, 1))
    threshold = kmeans.cluster_centers_.mean()
    labels = np.where(n_counts <= threshold, True, False)
    midsv_filtered_control = (m for m, label in zip(io.read_jsonl(path_midsv_control), labels) if label)
    io.write_jsonl(midsv_filtered_control, path_output)

    return path_output


def cache_normalized_indels(ARGS, path_midsv_sample: Path) -> None:
    allele, label, *_ = path_midsv_sample.stem.split("_")
    sequence = ARGS.fasta_alleles[allele]

    path_midsv_control = extract_path_control(ARGS, allele)
    path_midsv_n_filtered_control = extract_path_n_filtered_control(ARGS, path_midsv_control)

    cache_selected_control_by_similarity(
        ARGS, path_midsv_n_filtered_control, path_midsv_sample, path_midsv_sample.parent
    )

    path_midsv_similar_control = Path(path_midsv_sample.parent, f"{allele}_{label}_control.jsonl")

    _, indels_normalized_sample = summarize_indels(path_midsv_sample, sequence)
    _, indels_normalized_control = summarize_indels(path_midsv_similar_control, sequence)

    path_consensus = Path(ARGS.tempdir, ARGS.sample_name, "consensus")
    io.save_pickle(indels_normalized_sample, Path(path_consensus, f"{allele}_{label}_normalized_sample.pickle"))
    io.save_pickle(indels_normalized_control, Path(path_consensus, f"{allele}_{label}_normalized_control.pickle"))


def cache_mutation_loci(ARGS, clust_sample: list[dict]) -> None:
    # Separate clusters by label and cache them
    clust_sample.sort(key=lambda x: [x["ALLELE"], x["LABEL"]])
    path_consensus = Path(ARGS.tempdir, ARGS.sample_name, "consensus")
    for (allele, label), group in groupby(clust_sample, key=lambda x: [x["ALLELE"], x["LABEL"]]):
        io.write_jsonl(group, Path(path_consensus, f"{allele}_{label}_sample.jsonl"))

    # Cache normalized indels counts
    for path_midsv_sample in path_consensus.glob("*_sample.jsonl"):
        cache_normalized_indels(ARGS, path_midsv_sample)

    # Extract and cache mutation loci
    for path_indels_normalized_sample in path_consensus.glob("*_normalized_sample.pickle"):
        allele, label, *_ = path_indels_normalized_sample.stem.split("_")
        path_indels_normalized_control = Path(path_consensus, f"{allele}_{label}_normalized_control.pickle")
        sequence = ARGS.fasta_alleles[allele]
        path_knockin = Path(ARGS.tempdir, ARGS.sample_name, "knockin_loci", f"{allele}.pickle")

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

        io.save_pickle(mutation_loci, Path(path_consensus, f"{allele}_{label}_mutation_loci.pickle"))
