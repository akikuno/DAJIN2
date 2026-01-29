from __future__ import annotations

from itertools import groupby
from pathlib import Path

from DAJIN2.core.classification.allele_merger import merge_minor_alleles
from DAJIN2.utils import fileio


def calc_match(midsv_tags: str) -> int:
    """
    Calculate match score from MIDSV tags.
    1. Perfect match: 0
    2. Mismatch (insertion, deletion, inversion): -1 per event
    """
    mismatch_score = 0

    mismatch_score += midsv_tags.count("+")  # insertion
    mismatch_score += midsv_tags.count("-")  # deletion
    mismatch_score += midsv_tags.count("*")  # substitution
    mismatch_score += sum(tags.islower() for tags in midsv_tags)  # inversion

    return -mismatch_score


def score_allele(path_midsv: Path, allele: str) -> list[dict]:
    midsv_sample = fileio.read_jsonl(path_midsv)
    scored_alleles = []
    for dict_midsv in midsv_sample:
        score = calc_match(dict_midsv["MIDSV"])
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


def classify_alleles(TEMPDIR: Path, FASTA_ALLELES: dict, SAMPLE_NAME: str, no_filter: bool = False) -> list[dict]:
    score_of_each_alleles = []
    for allele in FASTA_ALLELES:
        path_midsv = Path(TEMPDIR, SAMPLE_NAME, "midsv", allele, f"{SAMPLE_NAME}_midsv.jsonl")
        score_of_each_alleles.extend(score_allele(path_midsv, allele))

    if no_filter:
        # Skip filtering when --no-filter is used
        score_of_each_alleles_merged = score_of_each_alleles
    else:
        # Apply existing filtering logic
        score_of_each_alleles_merged = merge_minor_alleles(score_of_each_alleles, threshold=5)

    return extract_alleles_with_max_score(score_of_each_alleles_merged)
