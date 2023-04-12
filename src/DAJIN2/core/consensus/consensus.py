from __future__ import annotations

import re
from itertools import groupby
from collections import defaultdict


def call_percentage(cssplits: list[str]) -> list[dict[str, float]]:
    """
    Call position weight matrix in defferent loci.
    Non defferent loci are annotated to "Match" or "Unknown(N)"
    """
    cssplits_transposed = [list(cs) for cs in zip(*cssplits)]
    coverage = len(cssplits)
    cons_percentage = []
    for cs_transposed in cssplits_transposed:
        count_cs = defaultdict(int)
        for cs in cs_transposed:
            count_cs[cs] += 1 / coverage * 100
        count_cs_sorted = dict(sorted(count_cs.items(), key=lambda x: x[1], reverse=True))
        cons_percentage.append(count_cs_sorted)
    return cons_percentage


def call_sequence(cons_percentage: list[dict[str, float]]) -> str:
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


def call_consensus(clust_sample: list[dict]) -> tuple[dict, dict]:
    cons_percentage = defaultdict(list)
    cons_sequence = defaultdict(list)
    clust_sample.sort(key=lambda x: x["LABEL"])
    for _, group in groupby(clust_sample, key=lambda x: x["LABEL"]):
        clust = list(group)
        keys = (
            clust[0]["ALLELE"],
            clust[0]["LABEL"],
            clust[0]["PERCENT"],
        )
        cssplits = [cs["CSSPLIT"].split(",") for cs in clust]
        cons_per = call_percentage(cssplits)
        cons_seq = call_sequence(cons_per)
        cons_percentage[keys] = cons_per
        cons_sequence[keys] = cons_seq
    return cons_percentage, cons_sequence


def detect_sv(cons_sequence: dict, threshold: int = 50) -> list[bool]:
    exists_sv = []
    for seq in cons_sequence.values():
        if "N" * threshold in seq:
            exists_sv.append(True)
        elif re.search(rf"(\+[ACGTN]\|){{{threshold}}}", seq):
            exists_sv.append(True)
        elif re.search(rf"(\-[ACGTN]){{{threshold}}}", seq):
            exists_sv.append(True)
        elif re.search(rf"(\*[ACGTN][ACGTN]){{{threshold}}}", seq):
            exists_sv.append(True)
        elif re.search(r"[acgtn]", seq):
            exists_sv.append(True)
        else:
            exists_sv.append(False)
    return exists_sv


def call_allele_name(cons_sequence: dict, FASTA_ALLELES: dict) -> dict[int, str]:
    exists_sv = detect_sv(cons_sequence)
    label_digits = len(str(len(cons_sequence)))
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

