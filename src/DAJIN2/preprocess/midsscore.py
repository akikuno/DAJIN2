"""
input: 各アレルごとのMIDS
    .tmpDAJIN/midsconv/barcode25_control.csv
    .tmpDAJIN/midsconv/barcode25_inversion.csv
    .tmpDAJIN/midsconv/barcode25_target.csv

1: .tmpDAJIN/midsconv/barcode25_{allele}.csvについて、各リードの先頭にアレル名を加える
2: .tmpDAJIN/midsconv/barcode25_{allele}.csvについて、スコアをつける
3: .tmpDAJIN/midsconv/barcode25_*について、もっともスコアが低いアレルのみを取り出す

output: 各リードごとに一番ミスマッチが少なかったアレル
    00001,target,100
    00002,control,221
    00003,control,81

例えば以下のリードがあった場合、
00001,target,"MDDDDM" → 2/6
00001,control,"M,M,M,M,D,S,3M" -> 5/10
00001,inversion,"M,4M,D,M" -> 2/
"""

import re
import pathlib
from itertools import groupby
from src.DAJIN2.utils.io import read_jsonl, write_jsonl


def calc(mids: str) -> float:
    mids_list = mids.split(",")
    mids_length = len(mids_list)
    for ins in mids_list:
        if ins[0].isdigit():
            mids_length += int(re.sub(r"M|D|S", "", ins))
    return mids_list.count("M") / mids_length


def extract_possible_allele_and_score(basename: str) -> list[dict]:
    readid_allele_score_of_all_alleles = []
    midspaths = pathlib.Path(".tmpDAJIN", "midsconv").glob(f"{basename}*")
    for midspath in midspaths:
        allele_name = midspath.stem.replace(f"{basename}_", "")
        for midscsv in midspath.read_text().split():
            readid = midscsv.split(",")[0]
            mids = ",".join(midscsv.split(",")[1:])
            score = calc(mids)
            readid_allele_score_of_all_alleles.append(
                {"readid": readid, "allele": allele_name, "score": score, "mids": mids}
            )
    readid_allele_score_of_all_alleles.sort(key=lambda x: x["readid"])
    readid_allele_score = []
    for _, group in groupby(readid_allele_score_of_all_alleles, key=lambda x: x["readid"]):
        max_score = 0
        for readinfo in group:
            if readinfo["score"] > max_score:
                max_read = readinfo
                max_score = readinfo["score"]
        readid_allele_score.append(max_read)
    return readid_allele_score


tmp = extract_possible_allele_and_score(sample_name)


#####
from collections import defaultdict

d = defaultdict(int)
for x in tmp:
    if x["allele"] == "control":
        print(x)
    d[x["allele"]] += 1

d
