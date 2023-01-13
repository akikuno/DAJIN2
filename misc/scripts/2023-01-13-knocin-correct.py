"""
課題
- Cables2 floxアレルにおいて、knock-in配列の補正がされていないため、knock-in配列上のシークエンスエラーによってアレルがクラスタリングされてしまっている
    - 1776番目の*GA変異
    - 2444番目の-T欠失
    - 2438の-C, 2456の-C, 2458の*TC
対策
- kcnok-in配列上の5merをとり、コントロールの5merと類似の配列を探す
    - knock-in配列を完全に含むように、周辺配列から1-2mer選択する
- 類似配列において、knock-in配列と変異プロファイルを比較し、変異の類似性が高い場合にはその変異を置き換える

結果
"""
import midsv
from pathlib import Path
from src.DAJIN2.core import preprocess, classification, clustering, consensus, report

# * flox insertion
SAMPLE, CONTROL, ALLELE, NAME, GENOME, DEBUG, THREADS = (
    "examples/flox-cables2/AyabeTask1/barcode31.fq.gz",
    "examples/flox-cables2/AyabeTask1/barcode42.fq.gz",
    "examples/flox-cables2/AyabeTask1/design_cables2.fa",
    "test-ayabe-task1",
    "mm10",
    True,
    14,
)


FASTA_ALLELES = preprocess.format_inputs.dictionize_allele(ALLELE)

allele = "flox"
sv = False
midsv_sample = midsv.read_jsonl(Path(TEMPDIR, "midsv", "barcode31_splice_flox.jsonl"))
midsv_control = midsv.read_jsonl(Path(TEMPDIR, "midsv", "barcode42_splice_flox.jsonl"))
midsv_sample[0]
midsv_control[3]
cssplits_sample = []
for samp in classif_sample:
    if samp["ALLELE"] == allele and not samp["SV"]:
        cssplits_sample.append(samp)

idx = 2444
idx = 1776
for cs in cssplits_sample:
    cs["CSSPLIT"].split(",")[idx]

cssplits_sample
sequence = ["TAGCTAAT"]
cssplits_sample = ["=T", "=A", "-A", "+G|C|T|=A", "A", "+C|-T", "=T"]

cssplits_decompose = []
for i, cssplits in enumerate(cssplits_sample):
    if cssplits.startswith("-") or cssplits == "N":
        continue


knockin_5mer = ["AGCTA", "TAACT"]

import difflib

difflib.get_close_matches(knockin_5mer, sequence, n=10, cutoff=0.0)

