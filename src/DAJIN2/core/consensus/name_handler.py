from __future__ import annotations

import re
from collections import defaultdict


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


def call_allele_name(
    cons_sequences: defaultdict[dict], cons_percentages: defaultdict[list], FASTA_ALLELES: dict, threshold: int = 50
) -> dict[int, str]:
    exists_sv = _detect_sv(cons_percentages, threshold)
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
