from importlib import reload
from pathlib import Path
from src.DAJIN2.core.preprocess import mappy_align

reload(mappy_align)


def test_to_sam_tyr():
    path_reference_fasta = Path("tests", "data", "mappy", "tyr_control.fa")
    path_query_fastx = Path("tests", "data", "mappy", "tyr_query.fq")
    value = mappy_align.to_sam(str(path_reference_fasta), str(path_query_fastx))
    value = list(value)
    answer = Path("tests", "data", "mappy", "tyr_query.sam").read_text().strip().split("\n")
    assert value == answer


def test_to_sam_stx2():
    path_reference_fasta = Path("tests", "data", "mappy", "stx2_control.fa")
    path_query_fastx = Path("tests", "data", "mappy", "stx2_del.fq")
    value = mappy_align.to_sam(str(path_reference_fasta), str(path_query_fastx))
    value = list(value)
    answer = Path("tests", "data", "mappy", "stx2_del.sam").read_text().strip().split("\n")
    assert value == answer


def test_to_sam_stx2_splice():
    path_reference_fasta = Path("tests", "data", "mappy", "stx2_control.fa")
    path_query_fastx = Path("tests", "data", "mappy", "stx2_splice.fq")
    value = mappy_align.to_sam(str(path_reference_fasta), str(path_query_fastx), preset="splice")
    value = list(value)
    answer = Path("tests", "data", "mappy", "stx2_del.sam").read_text().strip().split("\n")
    assert value == answer


def test_to_sam_threads():
    path_reference_fasta = Path("tests", "data", "mappy", "tyr_control.fa")
    path_query_fastx = Path("tests", "data", "mappy", "tyr_query.fq")
    value = mappy_align.to_sam(str(path_reference_fasta), str(path_query_fastx), THREADS=5)
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
