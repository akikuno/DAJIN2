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


def convert_to_percentage(cssplits: list[list[str]], mutation_loci: list[set[str]]) -> list[dict[str, float]]:
    """
    Convert sequences and mutations into percentages, annotating sequence errors as "SEQERROR".
    """
    # Transpose the cssplits to facilitate per-position processing
    cssplits_transposed = [list(cs) for cs in zip(*cssplits)]
    coverage = len(cssplits_transposed[0])

    cons_percentage = []
    for cs_transposed, mut_loci in zip(cssplits_transposed, mutation_loci):
        count_cs = defaultdict(float)
        for cs in cs_transposed:
            operator = cs[0]
            if operator in {"+", "-", "*"}:
                if operator not in mut_loci:
                    cs = "SEQERROR"
            count_cs[cs] += 1 / coverage * 100
        cons_percentage.append(dict(count_cs))

    return cons_percentage


def remove_all_n(cons_percentage: list[dict[str, float]]) -> list[dict[str, float]]:
    for c in cons_percentage:
        if len(c) == 1 and "N" in c:
            continue
        c.pop("N", None)

    return cons_percentage


def replace_sequence_error(cons_percentage: list[dict[str, float]]) -> list[dict[str, float]]:
    """
    Replace sequence error as distributing according to proportion of cs tags
    If a dictionary contains only "SEQERROR", it is replaced with {"N": 100}. Otherwise, "SEQERROR" is removed.
    """
    cons_percentage_replaced = []
    for cons_per in cons_percentage:
        # Replace a dictionary containing only "SEQERROR" with {"N": 100}
        if len(cons_per) == 1 and "SEQERROR" in cons_per:
            cons_percentage_replaced.append({"N": 100})
            continue
        cons_per.pop("SEQERROR", None)
        cons_percentage_replaced.append(cons_per)

    return cons_percentage_replaced


def adjust_to_100_percent(cons_percentage: list[dict[str, float]]) -> list[dict[str, float]]:
    adjusted_percentages = []

    for percentage_dict in cons_percentage:
        total = sum(percentage_dict.values())
        scaling_factor = 100 / total

        adjusted_dict = {key: value * scaling_factor for key, value in percentage_dict.items()}
        adjusted_percentages.append(adjusted_dict)

    return adjusted_percentages


def call_percentage(cssplits: list[list[str]], mutation_loci: list[set[str]]) -> list[dict[str, float]]:
    """call position weight matrix in defferent loci.
    - non defferent loci are annotated to "Match" or "Unknown(N)"
    """
    cons_percentage = convert_to_percentage(cssplits, mutation_loci)
    cons_percentage = remove_all_n(cons_percentage)
    cons_percentage = replace_sequence_error(cons_percentage)
    return adjust_to_100_percent(cons_percentage)


###########################################################
# main
###########################################################


@dataclass(frozen=True)
class ConsensusKey:
    allele: str
    label: int
    percent: float


def call_consensus(tempdir: Path, sample_name: str, clust_sample: list[dict]) -> tuple[dict[list], dict[str]]:
    clust_sample.sort(key=lambda x: [x["ALLELE"], x["LABEL"]])

    cons_percentages = {}
    cons_sequences = {}

    for (allele, label), group in groupby(clust_sample, key=lambda x: [x["ALLELE"], x["LABEL"]]):
        clust = list(group)

        path_consensus = Path(tempdir, sample_name, "consensus", allele, str(label))
        cons_mutation_loci = io.load_pickle(Path(path_consensus, "mutation_loci.pickle"))

        cssplits = [cs["CSSPLIT"].split(",") for cs in clust]
        cons_percentage = call_percentage(cssplits, cons_mutation_loci)

        key = ConsensusKey(allele, label, clust[0]["PERCENT"])
        cons_percentages[key] = cons_percentage
        cons_sequences[key] = call_sequence(cons_percentage)
    return cons_percentages, cons_sequences
