from collections.abc import Generator

import pytest

from src.DAJIN2.core.clustering.kmer_generator import generate_mutation_kmers


# Mock for fileio.read_jsonl
def mock_read_jsonl(path):
    return [{"MIDSV": "=A,-g,+t|t|t|=A,*ac,=N"}, {"MIDSV": "=A,-g,=T,*gt,=N"}]


@pytest.fixture
def mock_io_read_jsonl(monkeypatch):
    monkeypatch.setattr("DAJIN2.utils.fileio.read_jsonl", mock_read_jsonl)


# Test 1: Check if the function returns a generator
def test_is_generator(mock_io_read_jsonl):
    gen = generate_mutation_kmers("some_path", [set(), {"-"}, {"+", "*"}, {"*"}, set()])
    assert isinstance(gen, Generator)


# Test 2: Check the output with specific mutation_loci and compress_ins=True
def test_with_specific_mutation_loci_and_compress_ins_true(mock_io_read_jsonl):
    gen = generate_mutation_kmers("some_path", [set(), {"-"}, {"+", "*"}, {"*"}, set()], compress_ins=True)
    assert list(gen) == [
        ["@,@,@", "@,-g,+t|t|t|=A", "-g,+I=A,*ac", "+I=A,*ac,@", "@,@,@"],
        ["@,@,@", "@,-g,@", "@,@,@", "@,*gt,@", "@,@,@"],
    ]


# Test 3: Check the output with specific mutation_loci and compress_ins=False
def test_with_specific_mutation_loci_and_compress_ins_false(mock_io_read_jsonl):
    gen = generate_mutation_kmers("some_path", [set(), {"-"}, {"+", "*"}, {"*"}, set()], compress_ins=False)
    assert list(gen) == [
        ["@,@,@", "@,-g,+t|t|t|=A", "-g,+3=A,*ac", "+3=A,*ac,@", "@,@,@"],
        ["@,@,@", "@,-g,@", "@,@,@", "@,*gt,@", "@,@,@"],
    ]
