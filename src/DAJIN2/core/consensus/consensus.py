from __future__ import annotations

import pickle

from pathlib import Path
from typing import NamedTuple
from itertools import groupby
from collections import defaultdict

from DAJIN2.utils.cssplits_handler import find_n_boundaries

###########################################################
# call position weight matrix (cons_pergentage)
###########################################################


def remove_all_n(cons_percentage: list[dict[str, float]]) -> list[dict[str, float]]:
    for c in cons_percentage:
        if c == {"N": 100}:
            continue
        _ = c.pop("N", None)
    return cons_percentage


def replace_sequence_errror(cons_percentage: list[dict[str, float]]) -> list[dict[str, float]]:
    """replace sequence error as distributing according to proportion of cs tags"""
    cons_percentage_update = []
    for cons_per in cons_percentage:
        if "SEQERROR" not in cons_per:
            cons_percentage_update.append(cons_per)
            continue
        if len(cons_per) == 1 and cons_per["SEQERROR"]:
            cons_percentage_update.append({"N": 100})
            continue
        cons_per_update = dict()
        div = 100 / (sum(cons_per.values()) - cons_per["SEQERROR"])
        for key, val in cons_per.items():
            if key == "SEQERROR":
                continue
            cons_per_update[key] = val * div
        cons_percentage_update.append(cons_per_update)
    return cons_percentage_update


def adjust_to_100_percent(cons_percentage: list[dict[str, float]]) -> list[dict[str, float]]:
    adjusted_percentages = []

    for percentage_dict in cons_percentage:
        total = sum(percentage_dict.values())
        scaling_factor = 100 / total

        adjusted_dict = {key: value * scaling_factor for key, value in percentage_dict.items()}
        adjusted_percentages.append(adjusted_dict)

    return adjusted_percentages


def call_percentage(cssplits: list[str], mutation_loci) -> list[dict[str, float]]:
    """call position weight matrix in defferent loci.
    - non defferent loci are annotated to "Match" or "Unknown(N)"
    - sequence errors are annotated to "SEQERROR"
    """
    cssplits_transposed = [list(cs) for cs in zip(*cssplits)]
    coverage = len(cssplits)
    cons_percentage = []
    for cs_transposed, mut_loci in zip(cssplits_transposed, mutation_loci):
        count_cs = defaultdict(float)
        for cs in cs_transposed:
            if cs[0] in {"+", "-", "*"} and cs[0] not in mut_loci:
                cs = "SEQERROR"
            count_cs[cs] += 1 / coverage * 100
        cons_percentage.append(dict(count_cs))
    cons_percentage = remove_all_n(cons_percentage)
    cons_percentage = replace_sequence_errror(cons_percentage)
    return adjust_to_100_percent(cons_percentage)


###########################################################
# Call sequence
###########################################################


def _process_base(cons: str) -> str:
    if cons.startswith("="):  # match
        return cons.replace("=", "")
    elif cons.startswith("-"):  # deletion
        return ""
    elif cons.startswith("*"):  # substitution
        return cons[-1]
    elif cons.startswith("+"):  # insertion
        cons_ins = cons.split("|")
        if cons_ins[-1].startswith("="):  # match after insertion
            cons = cons.replace("=", "")
        elif cons_ins[-1].startswith("-"):  # deletion after insertion
            cons = "".join(cons_ins[:-1])
        elif cons_ins[-1].startswith("*"):  # substitution after insertion
            cons = "".join([*cons_ins[:-1], cons_ins[-1][-1]])
        return cons.replace("+", "").replace("|", "")
    return cons


def _call_sequence(cons_percentage: list[dict[str, float]]) -> str:
    consensus_sequence = []
    n_left, n_right = find_n_boundaries(cons_percentage)
    for i, cons_per in enumerate(cons_percentage):
        if n_left < i < n_right:
            cons = max(cons_per, key=cons_per.get)
            consensus_sequence.append(_process_base(cons))
        else:
            consensus_sequence.append("N")
    return "".join(consensus_sequence)


###########################################################
# main
###########################################################


class ConsensusKey(NamedTuple):
    allele: str
    label: int
    percent: float


def call_consensus(
    TEMPDIR: Path, SAMPLE_NAME: str, clust_sample: list[dict]
) -> tuple[defaultdict[list], defaultdict[str]]:
    cons_percentages = defaultdict(list)
    cons_sequences = defaultdict(str)
    clust_sample.sort(key=lambda x: x["LABEL"])
    for label, group in groupby(clust_sample, key=lambda x: x["LABEL"]):
        clust = list(group)
        allele = clust[0]["ALLELE"]
        key = ConsensusKey(allele, label, clust[0]["PERCENT"])
        cssplits = [cs["CSSPLIT"].split(",") for cs in clust]
        with open(Path(TEMPDIR, SAMPLE_NAME, "mutation_loci", f"{allele}.pickle"), "rb") as p:
            mutation_loci = pickle.load(p)
        cons_percentage = call_percentage(cssplits, mutation_loci)
        cons_percentages[key] = cons_percentage
        cons_sequences[key] = _call_sequence(cons_percentage)
    return cons_percentages, cons_sequences
