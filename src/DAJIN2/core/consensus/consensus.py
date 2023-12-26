from __future__ import annotations

from pathlib import Path
from typing import NamedTuple
from itertools import groupby
from collections import defaultdict

from DAJIN2.utils import io
from DAJIN2.utils.cssplits_handler import find_n_boundaries


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


def convert_consecutive_indels_to_match(cons_percentage: list[dict[str, float]]) -> list[dict[str, float]]:
    """
    Converts consecutive insertion and deletion operations on the same characters
    into a matched operation. This function processes a sequence of operations
    represented as a list of dictionaries. It cancels out an insertion ("+")
    operation followed downstream by a corresponding deletion ("-") operation on
    the same character. The function then consolidates the remaining operations
    and returns the updated sequence.
    """
    cons_percentage_reversed = cons_percentage[::-1].copy()

    i = 0
    while i < len(cons_percentage_reversed):
        current_cssplit = cons_percentage_reversed[i]
        current_consensus = max(current_cssplit, key=current_cssplit.get)

        if not current_consensus.startswith("+"):
            i += 1
            continue

        insertions = [base.lstrip("+") for base in current_consensus.split("|")[:-1]][::-1]

        # Extract deletions
        deletions = []
        count_del_min = float("inf")

        for j in range(1, len(insertions) + 1):
            if i + j >= len(cons_percentage_reversed):
                break

            next_consensus = max(cons_percentage_reversed[i + j], key=cons_percentage_reversed[i + j].get)

            if not next_consensus.startswith("-"):
                break

            deletions.append(next_consensus.lstrip("-"))
            count_del_min = min(count_del_min, cons_percentage_reversed[i + j][next_consensus])

        if insertions != deletions:
            i += 1
            continue

        # Format deletions
        count_ins = current_cssplit[current_consensus]

        for k, insertion in enumerate(insertions, 1):
            current_cssplit = cons_percentage_reversed[i + k]
            current_count_match = current_cssplit.get(f"={insertion}", 0)

            if count_del_min - count_ins == 0:
                current_cssplit[f"={insertion}"] = current_count_match + count_ins
                del current_cssplit[f"-{insertion}"]
            elif count_del_min - count_ins > 0:
                current_cssplit[f"={insertion}"] = current_count_match + count_ins
                current_cssplit[f"-{insertion}"] = count_del_min - count_ins
            else:
                current_cssplit[f"={insertion}"] = current_count_match + count_del_min
                current_cssplit[f"-{insertion}"] -= count_del_min
                if current_cssplit[f"-{insertion}"] <= 0:
                    del current_cssplit[f"-{insertion}"]

        # Format insertion
        cs_last = current_consensus.split("|")[-1]
        if count_del_min - count_ins >= 0:
            cons_percentage_reversed[i][cs_last] = count_ins
            del cons_percentage_reversed[i][current_consensus]
        else:
            cons_percentage_reversed[i][cs_last] = count_del_min
            cons_percentage_reversed[i][current_consensus] = count_ins - count_del_min

        i += len(insertions) + 1

    return cons_percentage_reversed[::-1]


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
    cons_percentage = convert_consecutive_indels_to_match(cons_percentage)
    return adjust_to_100_percent(cons_percentage)


###########################################################
# Call sequence
###########################################################


def cstag_to_base(cons: str) -> str:
    if cons.startswith("="):  # match
        return cons.replace("=", "")
    if cons.startswith("-"):  # deletion
        return ""
    if cons.startswith("*"):  # substitution
        return cons[-1]
    if cons.startswith("+"):  # insertion
        cons_ins = cons.split("|")
        if cons_ins[-1].startswith("="):  # match after insertion
            cons = cons.replace("=", "")
        elif cons_ins[-1].startswith("-"):  # deletion after insertion
            cons = "".join(cons_ins[:-1])
        elif cons_ins[-1].startswith("*"):  # substitution after insertion
            cons = "".join([*cons_ins[:-1], cons_ins[-1][-1]])
        return cons.replace("+", "").replace("|", "")
    return cons


def call_sequence(cons_percentage: list[dict[str, float]]) -> str:
    consensus_sequence = []
    n_left, n_right = find_n_boundaries(cons_percentage)
    for i, cons_per in enumerate(cons_percentage):
        if n_left < i < n_right:
            cons = max(cons_per, key=cons_per.get)
            consensus_sequence.append(cstag_to_base(cons))
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


def call_consensus(tempdir: Path, sample_name: str, clust_sample: list[dict]) -> tuple[dict[list], dict[str]]:
    cons_percentages = defaultdict(list)
    cons_sequences = defaultdict(str)

    path_consensus = Path(tempdir, sample_name, "consensus")

    clust_sample.sort(key=lambda x: [x["ALLELE"], x["LABEL"]])
    for (allele, label), group in groupby(clust_sample, key=lambda x: [x["ALLELE"], x["LABEL"]]):
        clust = list(group)

        mutation_loci = io.load_pickle(Path(path_consensus, f"{allele}_{label}_mutation_loci.pickle"))

        cssplits = [cs["CSSPLIT"].split(",") for cs in clust]
        cons_percentage = call_percentage(cssplits, mutation_loci)

        key = ConsensusKey(allele, label, clust[0]["PERCENT"])
        cons_percentages[key] = cons_percentage
        cons_sequences[key] = call_sequence(cons_percentage)
    return dict(cons_percentages), dict(cons_sequences)
