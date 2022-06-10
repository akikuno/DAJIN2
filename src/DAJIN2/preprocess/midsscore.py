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

import pathlib


def append_name(sample_name: str) -> dict:
    filePath = pathlib.Path(".tmpDAJIN", "midsconv").glob(f"{sample_name}*")
    for f in filePath:
        allele_name = f.name.replace(f"{sample_name}_", "").replace(".csv", "")
        tmp = f.read_text().split()
        scores = [calc(tmp)


def calc(mids: str) -> float:
    mids_list = mids.split(",")[1:]
    return mids_list.count("M") / len(mids_list)

