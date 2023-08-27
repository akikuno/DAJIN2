from __future__ import annotations

import midsv
from pathlib import Path
from itertools import groupby


def _calc_match(CSSPLIT: str) -> float:
    match_score = CSSPLIT.count("=")
    match_score -= CSSPLIT.count("+")  # insertion
    match_score -= sum(cs.islower() for cs in CSSPLIT)  # inversion
    cssplit = CSSPLIT.split(",")

    return match_score / len(cssplit)


def _score_allele(TEMPDIR: Path, allele: str, SAMPLE_NAME: str) -> list[dict]:
    midsv_sample = midsv.read_jsonl(Path(TEMPDIR, SAMPLE_NAME, "midsv", f"{allele}.json"))
    scored_alleles = []

    for dict_midsv in midsv_sample:
        score = _calc_match(dict_midsv["CSSPLIT"])
        dict_midsv.update({"SCORE": score, "ALLELE": allele})
        scored_alleles.append(dict_midsv)

    return scored_alleles


def _extract_alleles_with_max_score(score_of_each_alleles: list[dict]) -> list[dict]:
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
        score_of_each_alleles.extend(_score_allele(TEMPDIR, allele, SAMPLE_NAME))

    return _extract_alleles_with_max_score(score_of_each_alleles)
