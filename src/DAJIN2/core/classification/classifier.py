from __future__ import annotations

from pathlib import Path
from itertools import groupby
from collections import defaultdict

from DAJIN2.utils import io


def calc_match(cssplit: str) -> float:
    match_score = cssplit.count("=")
    match_score -= cssplit.count("+")  # insertion
    match_score -= sum(cs.islower() for cs in cssplit)  # inversion

    return match_score / len(cssplit.split(","))


def score_allele(path_midsv: Path, allele: str) -> list[dict]:
    midsv_sample = io.read_jsonl(path_midsv)

    scored_alleles = []
    for dict_midsv in midsv_sample:
        score = calc_match(dict_midsv["CSSPLIT"])
        dict_midsv.update({"SCORE": score, "ALLELE": allele})
        scored_alleles.append(dict_midsv)

    return scored_alleles


def extract_alleles_with_max_score(score_of_each_alleles: list[dict]) -> list[dict]:
    alleles_with_max_score = []
    score_of_each_alleles.sort(key=lambda x: x["QNAME"])
    for _, group in groupby(score_of_each_alleles, key=lambda x: x["QNAME"]):
        max_read = max(group, key=lambda x: x["SCORE"])
        del max_read["SCORE"]
        alleles_with_max_score.append(max_read)
    return alleles_with_max_score


##########################################################
# merge minor alleles
##########################################################


def count_alleles(score_of_each_alleles: list[dict]) -> dict[str, int]:
    score_of_each_alleles.sort(key=lambda x: x["QNAME"])
    allele_counts = defaultdict(int)
    for _, group in groupby(score_of_each_alleles, key=lambda x: x["QNAME"]):
        allele_predicted = max(group, key=lambda x: x["SCORE"])["ALLELE"]
        allele_counts[allele_predicted] += 1
    return allele_counts


# マイナーアレル, メジャーアレルの同定


def extract_minor_alleles(
    count_allele: dict[str, int], major_alleles: set[str], threshold_readnumber: int = 10
) -> dict[str, int]:
    return {
        allele: value
        for allele, value in count_allele.items()
        if allele not in major_alleles and value < threshold_readnumber
    }


def update_major_allele(
    count_allele: dict[str, int], major_alleles: set[str], threshold_readnumber: int = 10
) -> set[str]:
    new_major_alleles = {
        allele
        for allele, value in count_allele.items()
        if allele not in major_alleles and value >= threshold_readnumber
    }
    return major_alleles | new_major_alleles


# score_of_each_allelesから、マイナーアレルに該当するリードを抽出する


def extract_minor_groups(score_of_each_alleles: list[dict], minor_alleles: dict[str, int]) -> list[dict]:
    score_on_minor_alleles = []
    score_of_each_alleles.sort(key=lambda x: [x["QNAME"]])
    for _, group in groupby(score_of_each_alleles, key=lambda x: x["QNAME"]):
        group = list(group)
        allele_predicted = max(group, key=lambda x: x["SCORE"])["ALLELE"]
        if allele_predicted in minor_alleles:
            score_on_minor_alleles.extend(group)
    return score_on_minor_alleles


def merge_minor_alleles(score_of_each_alleles, threshold_readnumber: int = 10) -> list[dict]:
    score_of_each_alleles.sort(key=lambda x: [x["QNAME"]])
    major_alleles = set()
    minor_alleles = dict()

    allele_counts = count_alleles(score_of_each_alleles)
    minor_alleles = extract_minor_alleles(allele_counts, major_alleles, threshold_readnumber)
    major_alleles = update_major_allele(allele_counts, major_alleles, threshold_readnumber)

    minor_alleles_sorted = sorted(minor_alleles.items(), key=lambda x: x[1])

    score_on_minor_alleles = extract_minor_groups(score_of_each_alleles, minor_alleles)

    minor_allele_record = set()
    for minor_allele, _ in minor_alleles_sorted:
        minor_allele_record.add(minor_allele)
        score_on_minor_alleles.sort(key=lambda x: x["QNAME"])
        for _, group in groupby(score_on_minor_alleles, key=lambda x: x["QNAME"]):
            group = list(group)
            allele_predicted = max(group, key=lambda x: x["SCORE"])["ALLELE"]
            if allele_predicted != minor_allele:
                continue
            for g in group:
                if g["ALLELE"] in minor_allele_record:
                    g["SCORE"] = float("-inf")
        allele_counts = count_alleles(score_on_minor_alleles)
        minor_alleles = extract_minor_alleles(allele_counts, major_alleles, threshold_readnumber)
        major_alleles = update_major_allele(allele_counts, major_alleles, threshold_readnumber)
        score_on_minor_alleles = extract_minor_groups(score_on_minor_alleles, minor_alleles)
    return score_of_each_alleles


##########################################################
# main
##########################################################


def classify_alleles(TEMPDIR: Path, FASTA_ALLELES: dict, SAMPLE_NAME: str) -> list[dict]:
    score_of_each_alleles = []
    for allele in FASTA_ALLELES:
        path_midsv = Path(TEMPDIR, SAMPLE_NAME, "midsv", f"{allele}.json")
        score_of_each_alleles.extend(score_allele(path_midsv, allele))

    score_of_each_alleles_merged = merge_minor_alleles(score_of_each_alleles)

    return extract_alleles_with_max_score(score_of_each_alleles_merged)
