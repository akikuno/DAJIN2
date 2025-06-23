from __future__ import annotations

from pathlib import Path

from DAJIN2.core.preprocess.alignment import mapping


def test_to_sam_tyr():
    path_reference_fasta = Path("tests", "data", "preprocess", "mappy", "tyr_control.fa")
    path_query_fastx = Path("tests", "data", "preprocess", "mappy", "tyr_query.fq")
    test = mapping.to_sam(str(path_reference_fasta), str(path_query_fastx))
    test = list(test)
    answer = Path("tests", "data", "preprocess", "mappy", "tyr_query.sam").read_text().strip().split("\n")
    assert test == answer


def test_to_sam_stx2():
    path_reference_fasta = Path("tests", "data", "preprocess", "mappy", "stx2_control.fa")
    path_query_fastx = Path("tests", "data", "preprocess", "mappy", "stx2_del.fq")
    test = mapping.to_sam(str(path_reference_fasta), str(path_query_fastx))
    test = list(test)
    answer = Path("tests", "data", "preprocess", "mappy", "stx2_del.sam").read_text().strip().split("\n")
    assert test == answer


def test_to_sam_stx2_splice():
    path_reference_fasta = Path("tests", "data", "preprocess", "mappy", "stx2_control.fa")
    path_query_fastx = Path("tests", "data", "preprocess", "mappy", "stx2_splice.fq")
    test = mapping.to_sam(str(path_reference_fasta), str(path_query_fastx), preset="splice")
    test = list(test)
    answer = eval(Path("tests", "data", "preprocess", "mappy", "answer_stx2_splice.txt").read_text())
    assert test == answer


def test_to_sam_threads():
    path_reference_fasta = Path("tests", "data", "preprocess", "mappy", "tyr_control.fa")
    path_query_fastx = Path("tests", "data", "preprocess", "mappy", "tyr_query.fq")
    test = mapping.to_sam(str(path_reference_fasta), str(path_query_fastx), threads=5)
    test = list(test)
    answer = Path("tests", "data", "preprocess", "mappy", "tyr_query.sam").read_text().strip().split("\n")
    assert test == answer


########################################################################
# Create faidx
########################################################################


def test_make_faidx():
    path_fasta = "tests/data/preprocess/mappy/fasta.fa"
    test = mapping.make_faidx(path_fasta)
    answer = Path("tests/data/preprocess/mappy/fasta.fa.fai").read_text()
    assert test == answer


def test_make_faidx_wrap():
    path_fasta = "tests/data/preprocess/mappy/fasta_wrap.fa"
    test = mapping.make_faidx(path_fasta)
    answer = Path("tests/data/preprocess/mappy/fasta_wrap.fa.fai").read_text()
    assert test == answer


def test_make_faidx_real():
    path_fasta = "tests/data/preprocess/mappy/tyr_control.fa"
    test = mapping.make_faidx(path_fasta)
    answer = Path("tests/data/preprocess/mappy/tyr_control.fa.fai").read_text()
    assert test == answer
