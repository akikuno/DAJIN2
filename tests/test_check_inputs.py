import pytest
from src.DAJIN2.preprocess import check_inputs
from importlib import reload

reload(check_inputs)


def test_fastq_extension():
    with pytest.raises(AttributeError) as excinfo:
        check_inputs.fastq_extension("test.fqq")
    assert str(excinfo.value) == "test.fqq requires extensions either 'fastq', 'fastq.gz', 'fq' or 'fq.gz'"


# Contents check

fastq_path = "tests/data/empty.txt"


def test_fastq_error():
    with pytest.raises(AttributeError):
        _ = check_inputs.fastq_content(fastq_path)


def test_fasta_error():
    with pytest.raises(AttributeError):
        fasta_path = "tests/data/empty.txt"
        _ = check_inputs.fasta_content(fasta_path)


def test_fastq_without_error():
    fasta_path = "examples/nanosim/del-stx2/control.fq.gz"
    assert check_inputs.fastq_content(fasta_path) is None


def test_fasta_without_error():
    fasta_path = "examples/nanosim/del-stx2/design_stx2.fa"
    assert check_inputs.fasta_content(fasta_path) is None


def test_available_url_pass():
    url, flag_fail = check_inputs.available_url(["https://example.com"])
    assert ("https://example.com", False) == (url, flag_fail)


def test_available_url_fail():
    url, flag_fail = check_inputs.available_url(["https://example_xxx.com"])
    assert ("https://example_xxx.com", True) == (url, flag_fail)


def test_available_genome_pass():
    genome = "mm10"
    url = f"https://genome.ucsc.edu/cgi-bin/das/{genome}/dna?segment=1:1,10"
    assert check_inputs.available_genome("mm10", url) is None


def test_available_genome_fail():
    genome = "xxxx"
    url = f"https://genome.ucsc.edu/cgi-bin/das/{genome}/dna?segment=1:1,10"
    with pytest.raises(AttributeError) as excinfo:
        check_inputs.available_genome(genome, url)
    assert str(excinfo.value) == f"{genome} is not listed in UCSC genome browser"
