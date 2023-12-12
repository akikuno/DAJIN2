from __future__ import annotations

from pathlib import Path
from itertools import groupby

from DAJIN2.utils import io
from DAJIN2.core.preprocess.mutation_extractor import summarize_indels, extract_mutation_loci


def subset_clust(clust_sample: list[dict], num: int = 1000) -> list[dict]:
    clust_subset_sample = []
    clust_sample.sort(key=lambda x: x["LABEL"])
    for _, group in groupby(clust_sample, key=lambda x: x["LABEL"]):
        clust_subset_sample.extend(list(group)[:num])
    return clust_subset_sample


def cache_mutation_loci(ARGS, clust_sample: list[dict]) -> None:
    tempdir, sample_name, control_name, fasta_alleles = (
        ARGS.tempdir,
        ARGS.sample_name,
        ARGS.control_name,
        ARGS.fasta_alleles,
    )

    # Separate clusters by label and cache them
    path_consensus = Path(tempdir, sample_name, "consensus")
    clust_sample.sort(key=lambda x: [x["ALLELE"], x["LABEL"]])
    for (allele, label), group in groupby(clust_sample, key=lambda x: [x["ALLELE"], x["LABEL"]]):
        io.write_jsonl(group, Path(path_consensus, f"clust_{allele}_{label}.jsonl"))

    # Cache normalized indels counts
    for path_clust in path_consensus.glob("clust_*.jsonl"):
        _, allele, label = path_clust.stem.split("_")
        sequence = fasta_alleles[allele]
        _, indels_normalized = summarize_indels(path_clust, sequence)
        io.save_pickle(indels_normalized, Path(path_consensus, f"clust_{allele}_{label}_normalized.pickle"))

    # Extract and cache mutation loci
    path_mutation_control = Path(tempdir, control_name, "mutation_loci")
    for path_indels_normalized_sample in path_consensus.glob("clust_*.pickle"):
        _, allele, label, _ = path_indels_normalized_sample.stem.split("_")

        sequence = fasta_alleles[allele]

        file_name = f"{allele}_{sample_name}_normalized.pickle"
        if not Path(path_mutation_control, file_name).exists():
            file_name = f"{allele}_normalized.pickle"
        path_indels_normalized_control = Path(path_mutation_control, file_name)

        path_knockin = Path(tempdir, sample_name, "knockin_loci", f"{allele}.pickle")
        mutation_loci = extract_mutation_loci(
            sequence, path_indels_normalized_sample, path_indels_normalized_control, path_knockin
        )
        io.save_pickle(mutation_loci, Path(path_consensus, f"clust_{allele}_{label}_mutation_loci.pickle"))
