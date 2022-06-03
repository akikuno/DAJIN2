import pytest
from src.DAJIN2.preprocess.format_input import *

# Extension check
def test_check_fastq_extension():
    try:
        _ = check_fastq_extension("test.fqq")
    except:
        pytest.fail("error")


def test_check_fastq_extension():
    try:
        _ = check_fastq_extension("test.fq")
    except InputFileError as e:
        pytest.fail(e)


# Contents check

fastq_path = "tests/data/empty.txt"


def test_check_fastq_error():
    with pytest.raises(Exception):
        _ = check_fastq_content(fastq_path)


def test_check_fasta_error():
    with pytest.raises(Exception):
        fasta_path = "tests/data/empty.txt"
        _ = check_fasta_content(fasta_path)


def test_check_fastq_without_error():
    fasta_path = "examples/nanosim/del-stx2/control.fq.gz"
    assert check_fastq_content(fasta_path) == None


def test_check_fasta_without_error():
    fasta_path = "examples/nanosim/del-stx2/design_stx2.fa"
    assert check_fasta_content(fasta_path) == None

