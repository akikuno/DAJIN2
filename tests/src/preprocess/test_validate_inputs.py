import pytest
from DAJIN2.core.preprocess import validate_inputs
from importlib import reload

reload(validate_inputs)


###############################################################################
# validate File existance
###############################################################################


def test_exists():
    with pytest.raises(FileNotFoundError) as e:
        test = "filenotfound.txt"
        validate_inputs.exists(test)
    assert str(e.value) == f"{test} is not found"


###############################################################################
# validate FASTQ
###############################################################################


def test_fastq_extension():
    with pytest.raises(AttributeError) as e:
        test = "test.fqq"
        validate_inputs.fastq_extension("test.fqq")
    assert str(e.value) == f"{test} requires extensions either 'fastq', 'fastq.gz', 'fq' or 'fq.gz'"


def test_fastq_error_not_fastq_format():
    with pytest.raises(AttributeError):
        fastq_path = "tests/data/preprocess/validate_inputs/empty.fq"
        _ = validate_inputs.fastq_content(fastq_path)


def test_fastq_without_error():
    fasta_path = "examples/nanosim/del-stx2/control.fq.gz"
    assert validate_inputs.fastq_content(fasta_path) is None


###############################################################################
# validate FASTA
###############################################################################


def test_fasta_error_not_fasta_format():
    with pytest.raises(AttributeError) as e:
        fasta_path = "tests/data/preprocess/validate_inputs/empty.fa"
        _ = validate_inputs.fasta_content(fasta_path)
    assert str(e.value) == f"{fasta_path} is not a FASTA format"


def test_fasta_error_duplicated_identifiers():
    with pytest.raises(AttributeError) as e:
        fasta_path = "tests/data/preprocess/validate_inputs/duplicated_name.fa"
        _ = validate_inputs.fasta_content(fasta_path)
    assert str(e.value) == f"{fasta_path} must include unique identifiers"


def test_fasta_error_duplicated_sequences():
    with pytest.raises(AttributeError) as e:
        fasta_path = "tests/data/preprocess/validate_inputs/duplicated_seq.fa"
        _ = validate_inputs.fasta_content(fasta_path)
    assert str(e.value) == f"{fasta_path} must include unique DNA sequences"


def test_fasta_error_without_control():
    with pytest.raises(AttributeError) as e:
        fasta_path = "tests/data/preprocess/validate_inputs/no_control.fa"
        _ = validate_inputs.fasta_content(fasta_path)
    assert str(e.value) == f"{fasta_path} must include a 'control' sequence"


def test_fasta_without_error():
    fasta_path = "examples/nanosim/del-stx2/design_stx2.fa"
    assert validate_inputs.fasta_content(fasta_path) is None


###############################################################################
# validate URL
###############################################################################


@pytest.mark.skip("This test takes long time due to URL access")
def test_available_url_pass():
    flag = validate_inputs._check_url_availabilities(["https://example.com"])
    assert flag == [True]


@pytest.mark.skip("This test takes long time due to URL access")
def test_available_url_fail():
    flag = validate_inputs._check_url_availabilities(["https://example_xxx.com"])
    assert flag == [False]


@pytest.mark.skip("This test takes long time due to URL access")
def test_available_genome_pass():
    genome = "mm10"
    ucsc_url = "https://genome.ucsc.edu/"
    assert validate_inputs._is_listed(genome, ucsc_url) is None


@pytest.mark.skip("This test takes long time due to URL access")
def test_available_genome_fail():
    genome = "xxxx"
    ucsc_url = "https://genome.ucsc.edu/"
    with pytest.raises(AttributeError) as e:
        validate_inputs._is_listed(genome, ucsc_url)
    assert (
        str(e.value)
        == f"{genome} is not listed in UCSC genome browser. Available genomes are in {ucsc_url}/cgi-bin/das/dsn"
    )
