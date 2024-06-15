from __future__ import annotations

from collections import defaultdict
from itertools import groupby

##########################################################
# merge minor alleles
##########################################################


def count_allele_with_max_score(score_of_each_alleles: list[dict]) -> dict[str, int]:
    score_of_each_alleles.sort(key=lambda x: x["QNAME"])
    allele_counts = defaultdict(int)
    for _, group in groupby(score_of_each_alleles, key=lambda x: x["QNAME"]):
        allele_predicted = max(group, key=lambda x: x["SCORE"])["ALLELE"]
        allele_counts[allele_predicted] += 1

    return allele_counts


def extract_minor_alleles(allele_counts: dict[str, int], threshold_readnumber: int = 5) -> dict[str, int]:
    return {allele: value for allele, value in allele_counts.items() if value < threshold_readnumber}


def extract_major_alleles(allele_counts: dict[str, int], threshold_readnumber: int = 5) -> set[str]:
    return {allele for allele, value in allele_counts.items() if value >= threshold_readnumber}


def split_major_minor_alleles(
    score_of_each_alleles: list[dict], minor_alleles: dict[str, int]
) -> tuple[list[dict], list[dict]]:
    score_on_minor_alleles = []
    score_on_major_alleles = []
    score_of_each_alleles.sort(key=lambda x: [x["QNAME"]])
    for _, group in groupby(score_of_each_alleles, key=lambda x: x["QNAME"]):
        group = list(group)
        allele_predicted = max(group, key=lambda x: x["SCORE"])["ALLELE"]
        if allele_predicted in minor_alleles:
            score_on_minor_alleles.extend(group)
        else:
            score_on_major_alleles.extend(group)

    return score_on_major_alleles, score_on_minor_alleles


def extract_most_major_allele(major_alleles: dict[str, int]) -> str:
    return max(major_alleles, key=lambda x: major_alleles[x])


def replace_negative_inf_with_most_major_allele(
    score_of_minor_alleles: list[dict], all_allele_counts: dict[str, int]
) -> list[dict]:
    most_major_allele = extract_most_major_allele(all_allele_counts)
    score_of_minor_alleles = []
    for s in score_of_minor_alleles:
        if s["SCORE"] == float("-inf"):
            s["ALLELE"] = most_major_allele
        score_of_minor_alleles.append(s)
    return score_of_minor_alleles


def merge_minor_alleles(score_of_each_alleles, threshold_readnumber: int = 5) -> list[dict]:
    score_of_each_alleles.sort(key=lambda x: [x["QNAME"]])

    all_allele_counts = count_allele_with_max_score(score_of_each_alleles)
    minor_allele_counts = extract_minor_alleles(all_allele_counts, threshold_readnumber)

    score_of_alleles_merged, score_of_minor_alleles = split_major_minor_alleles(
        score_of_each_alleles, minor_allele_counts
    )

    score_of_minor_alleles.sort(key=lambda x: x["QNAME"])
    new_major_allele = set()
    for minor_allele in sorted(minor_allele_counts, key=minor_allele_counts.get):
        if minor_allele in new_major_allele:
            continue
        for _, group in groupby(score_of_minor_alleles, key=lambda x: x["QNAME"]):
            group = list(group)
            allele_predicted = max(group, key=lambda x: x["SCORE"])["ALLELE"]
            if allele_predicted != minor_allele:
                continue
            for g in group:
                if g["ALLELE"] in minor_allele:
                    g["SCORE"] = float("-inf")

        allele_counts = count_allele_with_max_score(score_of_minor_alleles)
        new_major_allele |= extract_major_alleles(allele_counts, threshold_readnumber)

    score_of_minor_alleles = replace_negative_inf_with_most_major_allele(score_of_minor_alleles, all_allele_counts)

    return score_of_alleles_merged + score_of_minor_alleles
