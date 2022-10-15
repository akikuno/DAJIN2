from importlib import reload
from pathlib import Path
from src.DAJIN2.core.preprocess import mappy_align

reload(mappy_align)


def test_to_sam():
    reffa = Path("tests", "data", "mappy", "tyr_control.fa")
    quefq = Path("tests", "data", "mappy", "tyr_query.fq")
    value = mappy_align.to_sam(str(reffa), str(quefq))
    value = list(value)
    answer = Path("tests", "data", "mappy", "tyr_query.sam").read_text().strip().split("\n")
    assert value == answer


########################################################################
# Create faidx
########################################################################


def test_make_faidx():
    path_fasta = "tests/data/mappy/fasta.fa"
    test = mappy_align.make_faidx(path_fasta)
    answer = Path("tests/data/mappy/fasta.fa.fai").read_text()
    assert test == answer


def test_make_faidx_wrap():
    path_fasta = "tests/data/mappy/fasta_wrap.fa"
    test = mappy_align.make_faidx(path_fasta)
    answer = Path("tests/data/mappy/fasta_wrap.fa.fai").read_text()
    assert test == answer


def test_make_faidx_real():
    path_fasta = "tests/data/mappy/tyr_control.fa"
    test = mappy_align.make_faidx(path_fasta)
    answer = Path("tests/data/mappy/tyr_control.fa.fai").read_text()
    assert test == answer
