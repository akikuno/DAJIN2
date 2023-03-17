"""
問題: 補正が甘いところがある
    - Albinoのアレルにおいて、2271番目のTの塩基が補正されていない
        - これは「T」TTTTTTと7塩基連続でTが続くので、`repetitive_del.py`か`sequence_error.py`のどちらかで補正されるべき
        - 2531番目のAの塩基も同じく補正されていない
対策
    - correct_repetitive_deletions.pyで"N"の塩基も補正の対象とした

問題: 上記の問題の後に、repetitiveではない箇所での変異において補正が甘いところがあることがわかった
    - 633の*GA および 1304 の-T など
原因
    - correct_sequence_error.pyにおいて、pvalueの計算のさいにリード数で正規化をしていなかった
    - そのため、リード数によってどの変異も有意差がついてしまっていた
対策
    - サンプルとコントロールのリード数で変異数を割ることで正規化し、pvalueを算出した
"""

from __future__ import annotations

import hashlib
from collections import defaultdict
from pathlib import Path
from importlib import reload
from src.DAJIN2.core import preprocess

import midsv

percent = 50
SAMPLE, CONTROL, ALLELE, NAME, GENOME, DEBUG, THREADS = (
    f"misc/data/tyr_albino_{percent}%.fq.gz",
    "misc/data/tyr_control.fq.gz",
    "misc/data/tyr_design.fasta",
    f"test-tyr-albino_{percent}%",
    "mm10",
    True,
    14,
)


midsv_sample = midsv.read_jsonl("DAJINResults/.tempdir/test-tyr-albino_50%/midsv/tyr_albino_50%_splice_albino.jsonl")

midsv_control = midsv.read_jsonl("DAJINResults/.tempdir/test-tyr-albino_50%/midsv/tyr_control_splice_albino.jsonl")

FASTA_ALLELES = preprocess.format_inputs.dictionize_allele(ALLELE)
allele = "albino"
sequence = FASTA_ALLELES[allele]
cssplits_sample = [m["CSSPLIT"].split(",") for m in midsv_sample]
cssplits_control = [m["CSSPLIT"].split(",") for m in midsv_control]

idx = 633
idx = 1304
count_sample = defaultdict(int)
for cs in cssplits_sample:
    count_sample[cs[idx]] += 1

count_control = defaultdict(int)
for cs in cssplits_control:
    count_control[cs[idx]] += 1

count_sample
count_control

