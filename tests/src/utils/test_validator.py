import pytest
from DAJIN2.utils import input_validator


###############################################################################
# validate File existance
###############################################################################


def test_exists():
    with pytest.raises(FileNotFoundError) as e:
        test = "filenotfound.txt"
        input_validator.validate_file_existence(test)
    assert str(e.value) == f"{test} is not found"


###############################################################################
# validate FASTQ
###############################################################################


def test_fastq_extension():
    with pytest.raises(ValueError) as e:
        test = "test.fqq"
        input_validator.validate_fastq_extension("test.fqq")
    assert str(e.value) == f"{test} requires extensions either 'fastq', 'fastq.gz', 'fq' or 'fq.gz'"


def test_fastq_error_not_fastq_format():
    with pytest.raises(ValueError):
        fastq_path = "tests/data/preprocess/validate_inputs/empty.fq"
        _ = input_validator.validate_fastq_content(fastq_path)


def test_fastq_without_error():
    fasta_path = "tests/data/preprocess/validate_inputs/control.fq.gz"
    assert input_validator.validate_fastq_content(fasta_path) is None


###############################################################################
# validate FASTA
###############################################################################


def test_non_proper_fasta_format():
    with pytest.raises(ValueError) as e:
        fasta_path = "tests/data/preprocess/validate_inputs/empty.fa"
        _ = input_validator.validate_fasta_content(fasta_path)
    assert str(e.value) == f"{fasta_path} is not a proper FASTA format"


def test_fasta_error_duplicated_identifiers():
    with pytest.raises(ValueError) as e:
        fasta_path = "tests/data/preprocess/validate_inputs/duplicated_name.fa"
        _ = input_validator.validate_fasta_content(fasta_path)
    assert str(e.value) == f"{fasta_path} must include unique identifiers"


def test_fasta_error_duplicated_sequences():
    with pytest.raises(ValueError) as e:
        fasta_path = "tests/data/preprocess/validate_inputs/duplicated_seq.fa"
        _ = input_validator.validate_fasta_content(fasta_path)
    assert str(e.value) == f"{fasta_path} must include unique DNA sequences"


def test_fasta_error_without_control():
    with pytest.raises(ValueError) as e:
        fasta_path = "tests/data/preprocess/validate_inputs/no_control.fa"
        _ = input_validator.validate_fasta_content(fasta_path)
    assert str(e.value) == f"One of the headers in the {fasta_path} must be '>control'"


def test_fasta_without_error():
    fasta_path = "tests/data/preprocess/validate_inputs/design_stx2.fa"
    assert input_validator.validate_fasta_content(fasta_path) is None


###############################################################################
# validate URL
###############################################################################

@pytest.mark.slow
def test_available_genome_pass():
    genome = "mm10"
    url_das = "https://genome.ucsc.edu/cgi-bin/das/dsn"
    assert input_validator.is_genome_in_ucsc_ids(genome, url_das) is True


@pytest.mark.slow
def test_available_genome_fail():
    genome = "mm12345"
    url_das = "https://genome.ucsc.edu/cgi-bin/das/dsn"
    assert input_validator.is_genome_in_ucsc_ids(genome, url_das) is False
