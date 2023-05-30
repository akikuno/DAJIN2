from __future__ import annotations

import re
from collections import defaultdict
from itertools import groupby


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
        count_cs_sorted = dict(sorted(count_cs.items(), key=lambda x: x[1], reverse=True))
        cons_percentage.append(count_cs_sorted)
    return cons_percentage


def _replace_percentage(cons_percentage: list[dict[str, float]]) -> list[dict[str, float]]:
    """replace sequence error as distributing according to proportion of cs tags"""
    cons_percentage_updated = []
    for dict_percentage in cons_percentage:
        if "SEQERROR" not in dict_percentage:
            cons_percentage_updated.append(dict_percentage)
            continue
        if len(dict_percentage) == 1 and dict_percentage["SEQERROR"]:
            cons_percentage_updated.append({"N": 100})
            continue
        dict_percentage_update = dict()
        div = 100 / (100 - dict_percentage["SEQERROR"])
        for key, val in dict_percentage.items():
            if key == "SEQERROR":
                continue
            dict_percentage_update[key] = round(val * div, 5)
        cons_percentage_updated.append(dict_percentage_update)
    return cons_percentage_updated


def _call_sequence(cons_percentage: list[dict[str, float]]) -> str:
    consensus_sequence = []
    for cons_per in cons_percentage:
        cons = max(cons_per, key=cons_per.get)
        if cons.startswith("="):
            cons = cons.replace("=", "")
        elif cons.startswith("-"):
            continue
        elif cons.startswith("*"):
            cons = cons[-1]
        elif cons.startswith("+"):
            cons_ins = cons.split("|")
            if cons_ins[-1].startswith("="):
                cons = cons.replace("=", "")
            elif cons_ins[-1].startswith("-"):
                cons = "".join(cons_ins[:-1])
            elif cons_ins[-1].startswith("*"):
                cons = "".join([*cons_ins[:-1], cons_ins[-1][-1]])
            cons = cons.replace("+", "")
            cons = cons.replace("|", "")
        consensus_sequence.append(cons)
    return "".join(consensus_sequence)


def _detect_sv(cons_percentage: defaultdict[list], threshold: int = 50) -> list[bool]:
    exists_sv = []
    for cons_per in cons_percentage.values():
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


def call_consensus(clust_sample: list[dict], MUTATION_LOCI_ALLELES) -> tuple[defaultdict[list], defaultdict[str]]:
    cons_percentage = defaultdict(list)
    cons_sequence = defaultdict(str)
    clust_sample.sort(key=lambda x: x["LABEL"])
    for _, group in groupby(clust_sample, key=lambda x: x["LABEL"]):
        clust = list(group)
        keys = (
            clust[0]["ALLELE"],
            clust[0]["LABEL"],
            clust[0]["PERCENT"],
        )
        mutation_loci = MUTATION_LOCI_ALLELES[clust[0]["ALLELE"]]
        cssplits = [cs["CSSPLIT"].split(",") for cs in clust]
        cons_per = _call_percentage(cssplits, mutation_loci)
        cons_per = _replace_percentage(cons_per)
        cons_seq = _call_sequence(cons_per)
        cons_percentage[keys] = cons_per
        cons_sequence[keys] = cons_seq
    return cons_percentage, cons_sequence


def call_allele_name(
    cons_sequence: defaultdict[dict], cons_percentage: defaultdict[list], FASTA_ALLELES: dict
) -> dict[int, str]:
    exists_sv = _detect_sv(cons_percentage)
    label_digits = len(str(len(cons_percentage)))
    allele_names = {}
    for is_sv, (keys, cons_seq) in zip(exists_sv, cons_sequence.items()):
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
