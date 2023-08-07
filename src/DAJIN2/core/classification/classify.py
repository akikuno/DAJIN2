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


###########################################################
# main
###########################################################


def classify_alleles(TEMPDIR: Path, FASTA_ALLELES: dict, SAMPLE_NAME: str) -> list[dict]:
    # Scoring
    score_of_each_alleles = []
    for allele in FASTA_ALLELES:
        midsv_sample = midsv.read_jsonl(Path(TEMPDIR, SAMPLE_NAME, "midsv", f"{allele}.json"))
        for dict_midsv in midsv_sample:
            score = _calc_match(dict_midsv["CSSPLIT"])
            dict_midsv.update({"SCORE": score})
            dict_midsv.update({"ALLELE": allele})
            score_of_each_alleles.append(dict_midsv)
    score_of_each_alleles.sort(key=lambda x: x["QNAME"])
    # Extract alleles with max scores
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
