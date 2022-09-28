import pytest
from src.DAJIN2.core.preprocess import format_inputs
from importlib import reload

reload(format_inputs)


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
