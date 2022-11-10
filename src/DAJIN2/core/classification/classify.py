from __future__ import annotations

from itertools import groupby
from pathlib import Path

import midsv


def calc_match(CSSPLIT: str) -> float:
    match_score = CSSPLIT.count("=")
    match_score -= CSSPLIT.count("+")  # insertion
    match_score -= sum(cs.islower() for cs in CSSPLIT)  # inversion
    return match_score / len(CSSPLIT.split(","))


def classify_alleles(paths_midsv: list[Path]) -> list[dict]:
    score_of_each_alleles = []
    for path_midsv in paths_midsv:
        allele = path_midsv.stem.split("_")[-1]
        for dict_midsv in midsv.read_jsonl(path_midsv):
            score = calc_match(dict_midsv["CSSPLIT"])
            dict_midsv.update({"SCORE": score})
            dict_midsv.update({"ALLELE": allele})
            score_of_each_alleles.append(dict_midsv)
    score_of_each_alleles.sort(key=lambda x: x["QNAME"])
    possible_allele = []
    for _, group in groupby(score_of_each_alleles, key=lambda x: x["QNAME"]):
        max_score = -float("inf")
        for readinfo in group:
            if readinfo["SCORE"] > max_score:
                max_score = readinfo["SCORE"]
                max_read = readinfo
                del max_read["SCORE"]
        possible_allele.append(max_read)
    return possible_allele
