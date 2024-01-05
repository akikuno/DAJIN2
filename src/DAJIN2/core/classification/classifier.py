from __future__ import annotations

from pathlib import Path
from itertools import groupby

from DAJIN2.utils import io
from DAJIN2.core.classification.allele_merger import merge_minor_alleles


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
# main
##########################################################


def classify_alleles(TEMPDIR: Path, FASTA_ALLELES: dict, SAMPLE_NAME: str) -> list[dict]:
    score_of_each_alleles = []
    for allele in FASTA_ALLELES:
        path_midsv = Path(TEMPDIR, SAMPLE_NAME, "midsv", f"{allele}.json")
        score_of_each_alleles.extend(score_allele(path_midsv, allele))

    score_of_each_alleles_merged = merge_minor_alleles(score_of_each_alleles)

    return extract_alleles_with_max_score(score_of_each_alleles_merged)
