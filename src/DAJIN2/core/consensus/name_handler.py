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


def _format_allele_label(label: int, total_labels: int) -> str:
    label_digits = len(str(total_labels))
    return f"{label:0{label_digits}}"


def _determine_suffix(cons_seq: str, fasta_allele: str, is_sv: bool) -> str:
    if cons_seq == fasta_allele:
        return "_intact"
    elif is_sv:
        return "_sv"
    else:
        return "_indels"


def _construct_allele_name(
    label: int, allele: str, cons_seq: str, fasta_allele: str, percent: int, is_sv: bool, total_labels: int
) -> str:
    label_format = _format_allele_label(label, total_labels)
    suffix = _determine_suffix(cons_seq, fasta_allele, is_sv)
    return f"allele{label_format}_{allele}{suffix}_{percent}%"


def call_allele_name(
    cons_sequences: dict[int, dict],
    cons_percentages: dict[int, list],
    FASTA_ALLELES: dict[str, str],
    threshold: int = 50,
) -> dict[int, str]:
    exists_sv = _detect_sv(cons_percentages, threshold)
    total_labels = len(cons_percentages)
    allele_names = {}

    for is_sv, (keys, cons_seq) in zip(exists_sv, cons_sequences.items()):
        allele_name = _construct_allele_name(
            keys.label, keys.allele, cons_seq, FASTA_ALLELES[keys.allele], keys.percent, is_sv, total_labels
        )
        allele_names[keys.label] = allele_name

    return allele_names


# def construct_allele_name(
#     cons_sequences: dict[NamedTuple, dict],
#     cons_percentages: dict[NamedTuple, list],
#     FASTA_ALLELES: dict[str, str],
#     threshold: int = 50,
# ) -> dict[int, str]:
#     exists_sv = _detect_sv(cons_percentages, threshold)
#     total_labels = len(cons_percentages)
#     allele_names = {}

#     for is_sv, (keys, cons_seq) in zip(exists_sv, cons_sequences.items()):
#         allele_name = construct_allele_name(
#             keys.label, keys.allele, cons_seq, FASTA_ALLELES[keys.allele], keys.percent, is_sv, total_labels
#         )
#         allele_names[keys.label] = allele_name

#     return allele_names


# def call_allele_name(
#     cons_sequences: defaultdict[dict], cons_percentages: defaultdict[list], FASTA_ALLELES: dict, threshold: int = 50
# ) -> dict[int, str]:
#     exists_sv = _detect_sv(cons_percentages, threshold)
#     label_digits = len(str(len(cons_percentages)))
#     allele_names = {}
#     for is_sv, (keys, cons_seq) in zip(exists_sv, cons_sequences.items()):
#         label_format = f"{keys.label:0{label_digits}}"
#         allele_name = f"allele{label_format}_{keys.allele}"
#         if cons_seq == FASTA_ALLELES[keys.allele]:
#             allele_name += "_intact"
#         elif is_sv:
#             allele_name += "_sv"
#         else:
#             allele_name += "_indels"
#         allele_name += f"_{keys.percent}%"
#         allele_names.update({keys.label: allele_name})
#     return allele_names


def update_key_by_allele_name(cons: dict, allele_names: dict[int, str]) -> dict:
    cons_update = dict()
    for key in cons:
        old_allele = cons[key]
        new_allele = allele_names[key]
        cons_update[new_allele] = old_allele
    return cons_update


def add_key_by_allele_name(clust_sample: list[dict], allele_names: dict[int, str]) -> list[dict]:
    for clust in clust_sample:
        label = clust["LABEL"]
        clust["NAME"] = allele_names[label]
    return clust_sample