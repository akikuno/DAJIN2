from __future__ import annotations

import midsv
from pathlib import Path
from itertools import groupby

# import re
# import json
# from collections import defaultdict
# from itertools import permutations
# from typing import Generator
# from DAJIN2.core.preprocess import align


# def _extract_diff_loci(TEMPDIR: Path) -> defaultdict[dict]:
#     """
#     Extract differencial loci between alleles
#         - The purpose is to lower match_score between very similar alleles such as point mutation.
#     """
#     fasta_alleles = list(Path(TEMPDIR, "fasta").iterdir())
#     fasta_alleles = [f for f in fasta_alleles if f.suffix != ".fai" or not re.search("insertion", f.suffix)]
#     mutation_alleles = defaultdict(dict)
#     for pair in list(permutations(fasta_alleles, 2)):
#         ref, query = pair
#         ref_allele = ref.stem
#         alignments = align.to_sam(ref, query, preset="splice")
#         alignments = list(alignments)
#         alignments = [a.split("\t") for a in alignments]
#         alignments_midsv = midsv.transform(alignments, midsv=False, cssplit=True, qscore=False)[0]
#         cssplits = alignments_midsv["CSSPLIT"].split(",")
#         mutations = dict()
#         mismatches = 0
#         for i, cs in enumerate(cssplits):
#             if cs.startswith("="):
#                 continue
#             mutations[i] = cs
#             mismatches += 1
#             mismatches += cs.count("|")  # insertion
#         if mismatches < 10:
#             mutation_alleles[ref_allele].update(mutations)
#     return mutation_alleles


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
