from __future__ import annotations

import pytest
from DAJIN2.utils import fastx_handler

########################################################################
# Extract basename
########################################################################


def test_extract_filename():
    fastq_path = "tests/hoge/fuga/example.fq.gz"
    test = fastx_handler.extract_filename(fastq_path)
    answer = "example"
    assert test == answer


def test_extract_filename_fail():
    fastq_path = "tests/hoge/fuga/.fq.gz"
    with pytest.raises(ValueError) as e:
        fastx_handler.extract_filename(fastq_path)
    assert str(e.value) == "Provided name is empty or consists only of whitespace"


def test_extract_filename_change_filename():
    fastq_path = "tests/hoge/fuga/x?y|z.fq.gz"
    test = fastx_handler.extract_filename(fastq_path)
    answer = "x-y-z"
    assert test == answer


########################################################################
# Convert allele file to dictionary type fasta format
########################################################################


def test_dictionize_allele_fasta():
    path_fasta = "tests/data/preprocess/format_input/fasta.fa"
    test = fastx_handler.dictionize_allele(path_fasta)
    answer = {"test1": "ACGTACGTACGTACGTACGTACGTACGTACGTACGT", "test2": "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"}
    assert test == answer


def test_dictionize_allele_fasta_wrap():
    path_fasta = "tests/data/preprocess/format_input/fasta_wrap.fa"
    test = fastx_handler.dictionize_allele(path_fasta)
    answer = {"test1": "ACGTACGTACGTACGTACGTACGTACGTACGTACGT", "test2": "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"}
    assert test == answer


def test_dictionize_allele_empty():
    path_fasta = "tests/data/preprocess/format_input/fasta_empty.fa"
    with pytest.raises(ValueError) as e:
        fastx_handler.dictionize_allele(path_fasta)
    assert str(e.value) == "Provided name is empty or consists only of whitespace"
