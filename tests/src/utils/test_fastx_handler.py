from __future__ import annotations

import pytest

from pathlib import Path


from DAJIN2.utils import fastx_handler


def test_sanitize_filename_with_valid_path():
    assert fastx_handler.sanitize_filename("valid_name") == "valid_name"
    assert fastx_handler.sanitize_filename(Path("valid/name")) == "valid-name"


def test_sanitize_filename_with_invalid_characters():
    assert fastx_handler.sanitize_filename("inva/lid:name?") == "inva-lid-name-"
    assert fastx_handler.sanitize_filename(Path("/inva>lid|name.")) == "-inva-lid-name-"


def test_sanitize_filename_with_whitespace():
    assert fastx_handler.sanitize_filename("  leading_space") == "leading_space"
    assert fastx_handler.sanitize_filename("trailing_space ") == "trailing_space-"


def test_sanitize_filename_with_empty_string():
    with pytest.raises(ValueError) as e:
        fastx_handler.sanitize_filename(" ")
    assert str(e.value) == "Provided FASTA/FASTQ is empty or consists only of whitespace"


def test_sanitize_filename_with_empty_path():
    with pytest.raises(ValueError) as e:
        fastx_handler.sanitize_filename("")
    assert str(e.value) == "Provided FASTA/FASTQ is empty or consists only of whitespace"


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
    assert str(e.value) == "Provided FASTA/FASTQ is empty or consists only of whitespace"


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
    assert str(e.value) == "Provided FASTA/FASTQ is empty or consists only of whitespace"
