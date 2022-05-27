import pytest
from src.DAJIN2.preprocess.format_input import *


fastq_path = "tests/data/empty.txt"


def test_check_fastq_error():
    with pytest.raises(Exception):
        _ = check_fastq(fastq_path)


def test_check_fasta_error():
    with pytest.raises(Exception):
        fasta_path = "tests/data/empty.txt"
        _ = check_fasta(fasta_path)


def test_check_fastq_without_error():
    fasta_path = "examples/nanosim/del-stx2/control.fq.gz"
    assert check_fastq(fasta_path) == None


def test_check_fasta_without_error():
    fasta_path = "examples/nanosim/del-stx2/design_stx2.fa"
    assert check_fasta(fasta_path) == None
