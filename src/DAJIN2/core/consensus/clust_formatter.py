from __future__ import annotations

from pathlib import Path
from itertools import groupby

from sklearn.cluster import MiniBatchKMeans

from DAJIN2.utils import io
from DAJIN2.core.preprocess.mutation_extractor import summarize_indels, extract_mutation_loci
from DAJIN2.core.consensus.similarity_searcher import cache_selected_control_by_similarity


def subset_clust(clust_sample: list[dict], num: int = 1000) -> list[dict]:
    clust_subset_sample = []
    clust_sample.sort(key=lambda x: x["LABEL"])
    for _, group in groupby(clust_sample, key=lambda x: x["LABEL"]):
        clust_subset_sample.extend(list(group)[:num])
    return clust_subset_sample


###########################################################
# cache_mutation_loci
###########################################################


def get_thresholds(path_indels_normalized_sample, path_indels_normalized_control) -> dict[str, float]:
    indels_normalized_sample = io.load_pickle(path_indels_normalized_sample)
    indels_normalized_control = io.load_pickle(path_indels_normalized_control)
    thresholds = dict()
    for mut in {"+", "-", "*"}:
        values_sample = indels_normalized_sample[mut]
        values_control = indels_normalized_control[mut]
        values_subtract = values_sample - values_control
        kmeans = MiniBatchKMeans(n_clusters=2, random_state=0).fit(values_subtract.reshape(-1, 1))
        threshold = kmeans.cluster_centers_.mean()
        thresholds[mut] = max(threshold, 0.05)
    return thresholds


def cache_normalized_indels(ARGS, path_midsv_sample: Path) -> None:
    allele, label, *_ = path_midsv_sample.stem.split("_")
    sequence = ARGS.fasta_alleles[allele]

    if Path(ARGS.tempdir, ARGS.control_name, "midsv", f"{allele}.json").exists():
        path_midsv_control = Path(ARGS.tempdir, ARGS.control_name, "midsv", f"{allele}.json")
    else:
        path_midsv_control = Path(ARGS.tempdir, ARGS.control_name, "midsv", f"{allele}_{ARGS.sample_name}.json")

    cache_selected_control_by_similarity(path_midsv_control, path_midsv_sample, path_midsv_sample.parent)

    path_midsv_filtered_control = Path(path_midsv_sample.parent, f"{allele}_{label}_control.jsonl")

    _, indels_normalized_sample = summarize_indels(path_midsv_sample, sequence)
    _, indels_normalized_control = summarize_indels(path_midsv_filtered_control, sequence)

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

        mutation_loci = extract_mutation_loci(
            sequence, path_indels_normalized_sample, path_indels_normalized_control, path_knockin, thresholds
        )

        io.save_pickle(mutation_loci, Path(path_consensus, f"{allele}_{label}_mutation_loci.pickle"))
