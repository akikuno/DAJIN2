from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from itertools import groupby
from pathlib import Path

from DAJIN2.utils import io
from DAJIN2.utils.cssplits_handler import call_sequence

###########################################################
# call position weight matrix (cons_pergentage)
###########################################################


def convert_to_percentage(
    cssplits: list[list[str]], mutation_loci: list[set[str]], sequence: str
) -> list[dict[str, float]]:
    """
    Convert sequences and mutations into percentages, annotating sequence errors as "SEQERROR".
    """
    # Transpose the cssplits to facilitate per-position processing
    cssplits_transposed = [list(cs) for cs in zip(*cssplits)]
    coverage = len(cssplits_transposed[0])

    cons_percentage = []
    for cs_transposed, mut_loci, nucleotide in zip(cssplits_transposed, mutation_loci, sequence):
        count_cs = defaultdict(float)
        for cs in cs_transposed:
            operator = cs[0]
            if operator in {"+", "-", "*"}:
                if operator not in mut_loci:
                    cs = f"={nucleotide.upper()}"
            count_cs[cs] += 1 / coverage * 100
        cons_percentage.append(dict(count_cs))

    return cons_percentage


def remove_all_n(cons_percentage: list[dict[str, float]]) -> list[dict[str, float]]:
    for c in cons_percentage:
        if len(c) == 1 and "N" in c:
            continue
        c.pop("N", None)

    return cons_percentage


def adjust_to_100_percent(cons_percentage: list[dict[str, float]]) -> list[dict[str, float]]:
    adjusted_percentages = []

    for percentage_dict in cons_percentage:
        total = sum(percentage_dict.values())
        scaling_factor = 100 / total

        adjusted_dict = {key: value * scaling_factor for key, value in percentage_dict.items()}
        adjusted_percentages.append(adjusted_dict)

    return adjusted_percentages


def call_percentage(cssplits: list[list[str]], mutation_loci: list[set[str]], sequence: str) -> list[dict[str, float]]:
    """Call position weight matrix in different loci. Non-different loci are annotated as "Match"."""
    cons_percentage = convert_to_percentage(cssplits, mutation_loci, sequence)
    cons_percentage = remove_all_n(cons_percentage)
    return adjust_to_100_percent(cons_percentage)


###########################################################
# main
###########################################################


@dataclass(frozen=True)
class ConsensusKey:
    allele: str
    label: int
    percent: float


def call_consensus(
    tempdir: Path, sample_name: str, fasta_alleles: dict[str, str], clust_sample: list[dict]
) -> tuple[dict[list], dict[str]]:
    clust_sample.sort(key=lambda x: [x["ALLELE"], x["LABEL"]])

    cons_percentages = {}
    cons_sequences = {}

    for (allele, label), group in groupby(clust_sample, key=lambda x: [x["ALLELE"], x["LABEL"]]):
        clust = list(group)
        sequence = fasta_alleles[allele]

        path_consensus = Path(tempdir, sample_name, "consensus", allele, str(label))
        cons_mutation_loci = io.load_pickle(Path(path_consensus, "mutation_loci.pickle"))

        cssplits = [cs["CSSPLIT"].split(",") for cs in clust]
        cons_percentage = call_percentage(cssplits, cons_mutation_loci, sequence)

        key = ConsensusKey(allele, label, clust[0]["PERCENT"])
        cons_percentages[key] = cons_percentage
        cons_sequences[key] = call_sequence(cons_percentage)
    return cons_percentages, cons_sequences
