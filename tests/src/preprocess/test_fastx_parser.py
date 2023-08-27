import pytest
from src.DAJIN2.core.preprocess import fastx_parser


########################################################################
# Extract basename
########################################################################


def test_extract_basename():
    fastq_path = "tests/hoge/fuga/example.fq.gz"
    test = fastx_parser.extract_basename(fastq_path)
    answer = "example"
    assert test == answer


def test_extract_basename_fail():
    fastq_path = "tests/hoge/fuga/.fq.gz"
    with pytest.raises(ValueError) as e:
        fastx_parser.extract_basename(fastq_path)
    assert str(e.value) == "Provided FASTA/FASTQ is empty or consists only of whitespace"


def test_extract_basename_change_filename():
    fastq_path = "tests/hoge/fuga/x?y|z.fq.gz"
    test = fastx_parser.extract_basename(fastq_path)
    answer = "x-y-z"
    assert test == answer


########################################################################
# Convert allele file to dictionary type fasta format
########################################################################


def test_dictionize_allele_fasta():
    path_fasta = "tests/data/preprocess/format_input/fasta.fa"
    test = fastx_parser.dictionize_allele(path_fasta)
    answer = {"test1": "ACGTACGTACGTACGTACGTACGTACGTACGTACGT", "test2": "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"}
    assert test == answer


def test_dictionize_allele_fasta_wrap():
    path_fasta = "tests/data/preprocess/format_input/fasta_wrap.fa"
    test = fastx_parser.dictionize_allele(path_fasta)
    answer = {"test1": "ACGTACGTACGTACGTACGTACGTACGTACGTACGT", "test2": "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"}
    assert test == answer


def test_dictionize_allele_empty():
    path_fasta = "tests/data/preprocess/format_input/fasta_empty.fa"
    with pytest.raises(ValueError) as e:
        fastx_parser.dictionize_allele(path_fasta)
    assert str(e.value) == "Provided FASTA/FASTQ is empty or consists only of whitespace"
