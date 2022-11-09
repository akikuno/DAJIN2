from __future__ import annotations
from collections import defaultdict, Counter
from pathlib import Path
from copy import deepcopy
from src.DAJIN2.core.clustering import find_knockin_loci
import midsv


def test_find_knockin_loci():
    TEMPDIR = Path("tests/data/find_knockin_loci")
    DICT_ALLELE = eval(Path("tests/data/find_knockin_loci/design_cables2.jsonl").read_text())
    CONTROL_NAME = "barcode42"
    knockin_loci = find_knockin_loci(TEMPDIR, DICT_ALLELE, CONTROL_NAME)
    test = knockin_loci["flox"]
    answer = eval(Path("tests/data/find_knockin_loci/answer.txt").read_text())
    assert test == answer
