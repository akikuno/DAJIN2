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


readid_allele_score = []

# def append_name(sample_name: str) -> dict:

midspaths = pathlib.Path(".tmpDAJIN", "midsconv").glob(f"{sample_name}*")
for midspath in midspaths:
    allele_name = midspath.stem.replace(f"{sample_name}_", "")
    for midscsv in midspath.read_text().split():
        readid = midscsv.split(",")[0]
        mids = ",".join(midscsv.split(",")[1:])
        score = calc(mids)
        readid_allele_score.append({"readid": readid, "allele": allele_name, "score": score, "mids": mids})

readid_allele_score.sort(key=lambda x: x["readid"])

readid_allele_score_extract = []
for key, group in groupby(readid_allele_score, key=lambda x: x["readid"]):
    max_score = 0
    for member in group:
        if member["score"] > max_score:
            tmp_read = member
            max_score = member["score"]
    readid_allele_score_extract.append(tmp_read)


#####
# from collections import defaultdict

# d = defaultdict(int)
# for x in readid_allele_score_extract:
#     d[x["allele"]] += 1

