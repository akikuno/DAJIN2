import pytest
from src.DAJIN2.core.preprocess import format_inputs
from importlib import reload

reload(format_inputs)

########################################################################
# Extract basename
########################################################################


def test_extract_basename():
    fastq_path = "tests/hoge/fuga/example.fq.gz"
    test = format_inputs.extract_basename(fastq_path)
    answer = "example"
    assert test == answer


def test_extract_basename_fail():
    fastq_path = "tests/hoge/fuga/.fq.gz"
    with pytest.raises(AttributeError) as e:
        format_inputs.extract_basename(fastq_path)
    assert str(e.value) == f"{fastq_path} is not valid file name"


def test_extract_basename_change_filename():
    fastq_path = "tests/hoge/fuga/x?y|z.fq.gz"
    test = format_inputs.extract_basename(fastq_path)
    answer = "x-y-z"
    assert test == answer


########################################################################
# Convert allele file to dictionary type fasta format
########################################################################


def test_dictionize_allele_fasta():
    path_fasta = "tests/data/format_input/fasta.fa"
    test = format_inputs.dictionize_allele(path_fasta)
    answer = {"test1": "ACGTACGTACGTACGTACGTACGTACGTACGTACGT", "test2": "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"}
    assert test == answer


def test_dictionize_allele_fasta_wrap():
    path_fasta = "tests/data/format_input/fasta_wrap.fa"
    test = format_inputs.dictionize_allele(path_fasta)
    answer = {"test1": "ACGTACGTACGTACGTACGTACGTACGTACGTACGT", "test2": "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"}
    assert test == answer


def test_dictionize_allele_empty():
    path_fasta = "tests/data/format_input/fasta_empty.fa"
    with pytest.raises(AttributeError) as e:
        format_inputs.dictionize_allele(path_fasta)
    assert str(e.value) == f"{path_fasta} contains an empty header"

