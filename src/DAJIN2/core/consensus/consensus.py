from __future__ import annotations

import re
import pickle
from pathlib import Path
from collections import defaultdict
from itertools import groupby

###########################################################
# call position weight matrix (cons_pergentage)
###########################################################


def _count_consecutive_n(cons_percentage: list[dict[str, float]]) -> tuple(int, int):
    seq = "".join(max(c, key=c.get) for c in cons_percentage)
    n_left = len(re.match(r"^N*", seq)[0]) if re.match(r"^N*", seq) else 0
    n_right = len(seq) - len(re.match(r"N*$", seq)[0]) if re.match(r"N*$", seq) else len(seq)
    return n_left, n_right


def _remove_nonconsecutive_n(cons_percentage: list[dict[str, float]]) -> list[dict[str, float]]:
    """remove nonconsecutive N"""
    cons_percentage_remove_n = []
    n_left, n_right = _count_consecutive_n(cons_percentage)
    for i, cons_per in enumerate(cons_percentage):
        if i < n_left or i >= n_right:
            cons_percentage_remove_n.append(cons_per)
            continue
        if not cons_per == {"N": 100}:
            cons_per.pop("N", None)
        cons_percentage_remove_n.append(cons_per)
    return cons_percentage_remove_n


def _update_percentage(cons_percentage: list[dict[str, float]]) -> list[dict[str, float]]:
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


def _call_percentage(cssplits: list[str], mutation_loci) -> list[dict[str, float]]:
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
        cons_percentage.append(count_cs)
    cons_percentage = _remove_nonconsecutive_n(cons_percentage)
    cons_percentage = _update_percentage(cons_percentage)
    return cons_percentage


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
    n_left, n_right = _count_consecutive_n(cons_percentage)
    for i, cons_per in enumerate(cons_percentage):
        if i < n_left or i >= n_right:
            consensus_sequence.append("N")
        else:
            cons = max(cons_per, key=cons_per.get)
            consensus_sequence.append(_process_base(cons))
    return "".join(consensus_sequence)


def _detect_sv(cons_percentages: defaultdict[list], threshold: int = 50) -> list[bool]:
    exists_sv = []
    for cons_per in cons_percentages.values():
        cons_cssplits = []
        for cssplit in cons_per:
            seq = max(cssplit, key=cssplit.get)
            cons_cssplits.append(seq)
        cons_cssplits = "".join(cons_cssplits)
        if "N" * threshold in cons_cssplits:
            exists_sv.append(True)
        elif re.search(rf"(\+[ACGTN]\|){{{threshold}}}", cons_cssplits):
            exists_sv.append(True)
        elif re.search(rf"(\-[ACGTN]){{{threshold}}}", cons_cssplits):
            exists_sv.append(True)
        elif re.search(rf"(\*[ACGTN][ACGTN]){{{threshold}}}", cons_cssplits):
            exists_sv.append(True)
        elif re.search(r"[acgtn]", cons_cssplits):
            exists_sv.append(True)
        else:
            exists_sv.append(False)
    return exists_sv


###########################################################
# main
###########################################################


def call_consensus(
    TEMPDIR: Path, SAMPLE_NAME: str, clust_sample: list[dict]
) -> tuple[defaultdict[list], defaultdict[str]]:
    cons_percentages = defaultdict(list)
    cons_sequences = defaultdict(str)
    clust_sample.sort(key=lambda x: x["LABEL"])
    for label, group in groupby(clust_sample, key=lambda x: x["LABEL"]):
        clust = list(group)
        allele = clust[0]["ALLELE"]
        keys = (allele, label, clust[0]["PERCENT"])
        with open(Path(TEMPDIR, SAMPLE_NAME, "mutation_loci", f"{allele}.pickle"), "rb") as p:
            mutation_loci = pickle.load(p)
        cssplits = [cs["CSSPLIT"].split(",") for cs in clust]
        cons_percentage = _call_percentage(cssplits, mutation_loci)
        # cons_percentage = _replace_percentage(cons_percentage)
        cons_percentages[keys] = cons_percentage
        cons_sequences[keys] = _call_sequence(cons_percentage)
    return cons_percentages, cons_sequences


# def call_consensus(clust_sample: list[dict], MUTATION_LOCI_LABELS) -> tuple[defaultdict[list], defaultdict[str]]:
#     cons_percentages = defaultdict(list)
#     cons_sequences = defaultdict(str)
#     clust_sample.sort(key=lambda x: x["LABEL"])
#     for label, group in groupby(clust_sample, key=lambda x: x["LABEL"]):
#         clust = list(group)
#         keys = (
#             clust[0]["ALLELE"],
#             clust[0]["LABEL"],
#             clust[0]["PERCENT"],
#         )
#         mutation_loci = MUTATION_LOCI_LABELS[label]
#         cssplits = [cs["CSSPLIT"].split(",") for cs in clust]
#         cons_percentage = _call_percentage(cssplits, mutation_loci)
#         cons_percentage = _replace_percentage(cons_percentage)
#         cons_percentages[keys] = cons_percentage
#         cons_sequences[keys] = _call_sequence(cons_percentage)
#     return cons_percentages, cons_sequences


def call_allele_name(
    cons_sequences: defaultdict[dict], cons_percentages: defaultdict[list], FASTA_ALLELES: dict
) -> dict[int, str]:
    exists_sv = _detect_sv(cons_percentages)
    label_digits = len(str(len(cons_percentages)))
    allele_names = {}
    for is_sv, (keys, cons_seq) in zip(exists_sv, cons_sequences.items()):
        ALLELE, LABEL, PERCENT = keys
        label_format = f"{LABEL:0{label_digits}}"
        allele_name = f"allele{label_format}_{ALLELE}"
        if cons_seq == FASTA_ALLELES[ALLELE]:
            allele_name += "_intact"
        elif is_sv:
            allele_name += "_sv"
        else:
            allele_name += "_indels"
        allele_name += f"_{PERCENT}%"
        allele_names.update({LABEL: allele_name})
    return allele_names


def update_key_by_allele_name(cons: dict, allele_names: dict[int, str]) -> dict:
    for key, allele_name in zip(list(cons.keys()), allele_names.values()):
        cons[allele_name] = cons.pop(key)
    return cons


def add_key_by_allele_name(clust_sample: list[dict], allele_names: dict[int, str]) -> list[dict]:
    for clust in clust_sample:
        label = clust["LABEL"]
        clust["NAME"] = allele_names[label]
    return clust_sample
