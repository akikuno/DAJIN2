import pytest
from src.DAJIN2.core.preprocess import check_inputs
from importlib import reload

reload(check_inputs)


###############################################################################
# Check File existance
###############################################################################


def test_exists():
    with pytest.raises(FileNotFoundError) as e:
        test = "filenotfound.txt"
        check_inputs.exists(test)
    assert str(e.value) == f"{test} is not found"


###############################################################################
# Check FASTQ
###############################################################################


def test_fastq_extension():
    with pytest.raises(AttributeError) as e:
        test = "test.fqq"
        check_inputs.fastq_extension("test.fqq")
    assert str(e.value) == f"{test} requires extensions either 'fastq', 'fastq.gz', 'fq' or 'fq.gz'"


def test_fastq_error_not_fastq_format():
    with pytest.raises(AttributeError):
        fastq_path = "tests/data/check_inputs/empty.fq"
        _ = check_inputs.fastq_content(fastq_path)


def test_fastq_without_error():
    fasta_path = "examples/nanosim/del-stx2/control.fq.gz"
    assert check_inputs.fastq_content(fasta_path) is None


###############################################################################
# Check FASTA
###############################################################################


def test_fasta_error_not_fasta_format():
    with pytest.raises(AttributeError) as e:
        fasta_path = "tests/data/check_inputs/empty.fa"
        _ = check_inputs.fasta_content(fasta_path)
    assert str(e.value) == f"{fasta_path} is not a FASTA format"


def test_fasta_error_duplicated_identifiers():
    with pytest.raises(AttributeError) as e:
        fasta_path = "tests/data/check_inputs/duplicated_name.fa"
        _ = check_inputs.fasta_content(fasta_path)
    assert str(e.value) == f"{fasta_path} must include unique identifiers"


def test_fasta_error_duplicated_sequences():
    with pytest.raises(AttributeError) as e:
        fasta_path = "tests/data/check_inputs/duplicated_seq.fa"
        _ = check_inputs.fasta_content(fasta_path)
    assert str(e.value) == f"{fasta_path} must include unique DNA sequences"


def test_fasta_error_without_control():
    with pytest.raises(AttributeError) as e:
        fasta_path = "tests/data/check_inputs/no_control.fa"
        _ = check_inputs.fasta_content(fasta_path)
    assert str(e.value) == f"{fasta_path} must include a 'control' sequence"


def test_fasta_without_error():
    fasta_path = "examples/nanosim/del-stx2/design_stx2.fa"
    assert check_inputs.fasta_content(fasta_path) is None


###############################################################################
# Check URL
###############################################################################


def test_available_url_pass():
    url, flag_fail = check_inputs.available_url(["https://example.com"])
    assert ("https://example.com", False) == (url, flag_fail)


def test_available_url_fail():
    url, flag_fail = check_inputs.available_url(["https://example_xxx.com"])
    assert ("https://example_xxx.com", True) == (url, flag_fail)


def test_available_genome_pass():
    genome = "mm10"
    ucsc_url = "https://genome.ucsc.edu/"
    assert check_inputs.available_genome(genome, ucsc_url) is None


def test_available_genome_fail():
    genome = "xxxx"
    ucsc_url = "https://genome.ucsc.edu/"
    with pytest.raises(AttributeError) as e:
        check_inputs.available_genome(genome, ucsc_url)
    assert (
        str(e.value)
        == f"{genome} is not listed in UCSC genome browser. Available genomes are in {ucsc_url}/cgi-bin/das/dsn"
    )
