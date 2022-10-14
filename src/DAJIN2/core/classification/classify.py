from __future__ import annotations
import midsv
from pathlib import Path
from itertools import groupby


def calc_match(CSSPLIT: str) -> float:
    match_score = CSSPLIT.count("=")
    match_score -= CSSPLIT.count("+")  # insertion
    match_score -= sum(cs.islower() for cs in CSSPLIT)  # inversion
    return match_score / len(CSSPLIT.split(","))


def classify_alleles(path_midsv: Path, sample_name: str) -> list[dict]:
    score_of_each_allels = []
    for midsvpath in path_midsv:
        ALLELE = midsvpath.stem.replace(f"{sample_name}_", "")
        for midsvcsv in midsv.read_jsonl(midsvpath):
            QNAME = midsvcsv["QNAME"]
            CSSPLIT = midsvcsv["CSSPLIT"]
            SCORE = calc_match(CSSPLIT)
            score_of_each_allels.append({"QNAME": QNAME, "ALLELE": ALLELE, "SCORE": SCORE, "CSSPLIT": CSSPLIT})
    score_of_each_allels.sort(key=lambda x: x["QNAME"])
    possible_allele = []
    for _, group in groupby(score_of_each_allels, key=lambda x: x["QNAME"]):
        max_score = -float("inf")
        for readinfo in group:
            if readinfo["SCORE"] > max_score:
                max_score = readinfo["SCORE"]
                max_read = readinfo
                del max_read["SCORE"]
        possible_allele.append(max_read)
    return possible_allele
