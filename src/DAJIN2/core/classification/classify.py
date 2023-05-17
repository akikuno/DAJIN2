from __future__ import annotations

from collections import defaultdict
from itertools import groupby, permutations
from pathlib import Path

import midsv
import json
from typing import Generator
from DAJIN2.core.preprocess import mappy_align


def _extract_diff_loci(TEMPDIR: Path) -> defaultdict[dict]:
    """
    Extract differencial loci between alleles
        - The purpose is to lower match_score between very similar alleles such as point mutation.
    """
    fasta_alleles = list(Path(TEMPDIR, "fasta").iterdir())
    fasta_alleles = [f for f in fasta_alleles if f.suffix != ".fai"]
    mutation_alleles = defaultdict(dict)
    for pair in list(permutations(fasta_alleles, 2)):
        ref, query = pair
        ref_allele = ref.stem
        alignments = mappy_align.to_sam(ref, query, preset="splice")
        alignments = list(alignments)
        alignments = [a.split("\t") for a in alignments]
        alignments_midsv = midsv.transform(alignments, midsv=False, cssplit=True, qscore=False)[0]
        cssplits = alignments_midsv["CSSPLIT"].split(",")
        mutations = dict()
        for i, cs in enumerate(cssplits):
            if cs.startswith("="):
                continue
            mutations.update({i: cs})
        if len(mutations) < 10:
            mutation_alleles[ref_allele].update(mutations)
    return mutation_alleles


def _calc_match(CSSPLIT: str, mutations: dict) -> float:
    match_score = CSSPLIT.count("=")
    match_score -= CSSPLIT.count("+")  # insertion
    match_score -= sum(cs.islower() for cs in CSSPLIT)  # inversion
    cssplit = CSSPLIT.split(",")
    for i, mut in mutations.items():
        if cssplit[i] == mut:
            match_score = 0
    return match_score / len(cssplit)


def read_json(filepath) -> Generator[dict[str, str]]:
    with open(filepath, "r") as f:
        for line in f:
            yield json.loads(line)


###########################################################
# main
###########################################################


def classify_alleles(TEMPDIR: Path, FASTA_ALLELES: dict, SAMPLE_NAME: str) -> list[dict]:
    mutations = _extract_diff_loci(TEMPDIR)
    # Scoring
    score_of_each_alleles = []
    for allele in FASTA_ALLELES:
        midsv_sample = read_json(Path(TEMPDIR, "midsv", f"{SAMPLE_NAME}_{allele}.json"))
        for dict_midsv in midsv_sample:
            score = _calc_match(dict_midsv["CSSPLIT"], mutations[allele])
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
