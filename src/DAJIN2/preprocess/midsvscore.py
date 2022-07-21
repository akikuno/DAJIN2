from __future__ import annotations
import midsv
from pathlib import Path
from itertools import groupby


def calc_match(CSSPLIT: str) -> float:
    match_score = CSSPLIT.count("=")
    match_score -= CSSPLIT.count("N")  # unknown
    match_score -= CSSPLIT.count("+")  # insertion
    match_score -= sum(cs.islower() for cs in CSSPLIT)  # inversion
    return match_score / len(CSSPLIT.split(","))


def extract_possible_allele_and_score(sample_name: str) -> list[dict]:
    score_of_each_allels = []
    midsvpaths = list(Path(".tmpDAJIN", "midsv").glob(f"{sample_name}*"))
    # allele_num = len(midsvpaths)
    for midsvpath in midsvpaths:
        allele_name = midsvpath.stem.replace(f"{sample_name}_", "")
        for midsvcsv in midsv.read_jsonl(midsvpath):
            QNAME = midsvcsv["QNAME"]
            CSSPLIT = midsvcsv["CSSPLIT"]
            score = calc_match(CSSPLIT)
            score_of_each_allels.append({"QNAME": QNAME, "ALLELE": allele_name, "SCORE": score, "CSSPLIT": CSSPLIT})
    score_of_each_allels.sort(key=lambda x: x["QNAME"])
    possible_allele_and_score = []
    for _, group in groupby(score_of_each_allels, key=lambda x: x["QNAME"]):
        group = list(group)
        # if allele_num != len(group):
        #     continue
        max_score = 0
        for readinfo in group:
            if readinfo["SCORE"] > max_score:
                max_read = readinfo
                max_score = readinfo["SCORE"]
        possible_allele_and_score.append(max_read)
    return possible_allele_and_score


#####
# tmp = extract_possible_allele_and_its_score(control_name)
# from collections import defaultdict

# d = defaultdict(int)
# for x in tmp:
#     if x["ALLELE"] == "inversion":
#         print(x["QNAME"], x["SCORE"], x["CSSPLIT"])
#     d[x["ALLELE"]] += 1
# d
