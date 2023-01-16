from __future__ import annotations

"""
課題
- Cables2 floxアレルにおいて、knock-in配列の補正がされていないため、knock-in配列上のシークエンスエラーによってアレルがクラスタリングされてしまっている
    - 1776番目の*GA変異
    - 2444番目の-T欠失
    - 2438の-C, 2456の-C, 2458の*TC
対策
- kcnok-in配列上の5merをとり、コントロールの5merと類似の配列を探す
    - FASTA_ALLELEのマッピングによりknock-in配列箇所を探す
    - knock-in配列を完全に含むように、周辺配列から1-2mer選択する
- 類似配列において、knock-in配列と変異プロファイルを比較し、変異の類似性が高い場合にはその変異を置き換える

メモ
- knock-in配列のスタート
    left-loxp: 1733
    right-loxp: 2428
結果
"""
import midsv
from collections import defaultdict
from pathlib import Path
from src.DAJIN2.core import preprocess, classification, clustering, consensus, report
from itertools import groupby
from itertools import permutations
from collections import defaultdict
from pathlib import Path

from src.DAJIN2.core.preprocess import mappy_align

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

TEMPDIR = Path("DAJINResults", ".tempdir", NAME)
FASTA_ALLELES = preprocess.format_inputs.dictionize_allele(ALLELE)

allele = "flox"
sv = False
midsv_sample = midsv.read_jsonl(Path(TEMPDIR, "midsv", f"barcode31_splice_{allele}.jsonl"))
midsv_control = midsv.read_jsonl(Path(TEMPDIR, "midsv", f"barcode42_splice_{allele}.jsonl"))
midsv_sample[1]
midsv_control[3]
# cssplits_sample = []
# for samp in classif_sample:
#     if samp["ALLELE"] == allele and not samp["SV"]:
#         cssplits_sample.append(samp)

idx = 2444
idx = 1776
# left_start = 1733
count_sample = defaultdict(int)
for cs in midsv_sample:
    x = cs["CSSPLIT"].split(",")[idx]
    count_sample[x] += 1

count_control = defaultdict(int)
for cs in midsv_control:
    x = cs["CSSPLIT"].split(",")[idx]
    count_control[x] += 1

count_sample
count_control

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


def extract_diff_loci(TEMPDIR) -> defaultdict[dict]:
    """
    Extract differencial loci between alleles
        - The purpose is to lower match_score between very similar alleles such as point mutation.
    """
    fasta_alleles = list(Path(TEMPDIR, "fasta").iterdir())
    mutation_alleles = defaultdict(dict)
    for comb in list(permutations(fasta_alleles, 2)):
        ref, query = comb
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

