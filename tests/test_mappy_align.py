from importlib import reload
from pathlib import Path
from src.DAJIN2.preprocess import mappy_align

reload(mappy_align)


def test_to_sam():
    reffa = Path("tests", "data", "mappy", "tyr_control.fa")
    quefq = Path("tests", "data", "mappy", "tyr_query.fq")
    value = mappy_align.to_sam(str(reffa), str(quefq))
    value = list(value)
    answer = Path("tests", "data", "mappy", "tyr_query.sam").read_text().strip().split("\n")
    assert value == answer
