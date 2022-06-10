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


def calc(mids: str) -> float:
    mids_list = mids.split(",")[1:]
    mids_length = len(mids_list)
    for ins in mids_list:
        if ins[0].isdigit():
            mids_length += int(re.sub(r"M|D|S", "", ins))
    return mids_list.count("M") / mids_length


readid_allele_score = []

# def append_name(sample_name: str) -> dict:

midspath = pathlib.Path(".tmpDAJIN", "midsconv").glob(f"{sample_name}*")
for mids in midspath:
    allele_name = mids.stem.replace(f"{sample_name}_", "")
    for midscsv in mids.read_text().split():
        readid = midscsv.split(",")[0]
        score = calc(midscsv)
        readid_allele_score.append({"readid": readid, "allele": allele_name, "score": score})

readid_allele_score.sort(key=lambda x: x["readid"])

readid_allele_score_extract = []
for key, group in groupby(readid_allele_score, key=lambda x: x["readid"]):
    max_score = 0
    for member in group:
        if member["score"] > max_score:
            tmp_read = member
            max_score = member["score"]
    readid_allele_score_extract.append(tmp_read)

from collections import defaultdict

d = defaultdict(int)
for x in readid_allele_score_extract:
    d[x["allele"]] += 1
