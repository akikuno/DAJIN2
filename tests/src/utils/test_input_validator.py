from __future__ import annotations

from pathlib import Path

import pytest
from DAJIN2.utils import input_validator

###############################################################################
# validate File existance and the extentions
###############################################################################


def test_exists():
    with pytest.raises(FileNotFoundError) as e:
        test = "filenotfound.txt"
        input_validator.validate_file_existence(test)
    assert str(e.value) == f"{test} is not found"


def test_return_file_extension():
    with pytest.raises(ValueError) as e:
        test = Path("test.fqq")
        expected = f"{test} requires extensions either .fastq, .fastq.gz, .fq, .fq.gz, .fasta, .fasta.gz, .fa, .fa.gz, or .bam"
        input_validator.return_file_extension(test)
    assert str(e.value) == expected


###############################################################################
# validate FASTQ
###############################################################################


def test_validate_fastq_content_empty():
    with pytest.raises(ValueError):
        fastq_path = "tests/data/utils/validate_inputs/empty.fq"
        _ = input_validator.validate_fastq_content(fastq_path)


def test_validate_fastq_content_without_error():
    fasta_path = "tests/data/utils/validate_inputs/control.fq.gz"
    assert input_validator.validate_fastq_content(fasta_path) is None


###############################################################################
# validate FASTA
###############################################################################


def test_non_proper_fasta_format():
    with pytest.raises(ValueError) as e:
        fasta_path = "tests/data/utils/validate_inputs/empty.fa"
        _ = input_validator.validate_fasta_content(fasta_path)
    assert str(e.value) == f"{fasta_path} is not a proper FASTA format"


def test_validate_fasta_content_no_seq():
    with pytest.raises(ValueError):
        fasta_path = "tests/data/utils/validate_inputs/no_seq.fa"
        _ = input_validator.validate_fasta_content(fasta_path)


def test_fasta_error_duplicated_identifiers():
    with pytest.raises(ValueError) as e:
        fasta_path = "tests/data/utils/validate_inputs/duplicated_name.fa"
        _ = input_validator.validate_fasta_content(fasta_path)
    assert str(e.value) == f"{fasta_path} must include unique identifiers"


def test_fasta_error_duplicated_sequences():
    with pytest.raises(ValueError) as e:
        fasta_path = "tests/data/utils/validate_inputs/duplicated_seq.fa"
        _ = input_validator.validate_fasta_content(fasta_path)
    assert str(e.value) == f"{fasta_path} must include unique DNA sequences"


def test_fasta_error_without_control():
    with pytest.raises(ValueError) as e:
        fasta_path = "tests/data/utils/validate_inputs/no_control.fa"
        _ = input_validator.validate_fasta_content(fasta_path, allele_file=True)
    assert str(e.value) == f"One of the headers in the {fasta_path} must be '>control'"


def test_fasta_without_error():
    fasta_path = "tests/data/utils/validate_inputs/design_stx2.fa"
    assert input_validator.validate_fasta_content(fasta_path) is None


###############################################################################
# validate URL
###############################################################################

server_lists = {
    "blat": [
        "https://genome.ucsc.edu/cgi-bin/hgBlat",
        "https://genome-asia.ucsc.edu/cgi-bin/hgBlat",
        "https://genome-euro.ucsc.edu/cgi-bin/hgBlat",
    ],
    "das": [
        "https://genome.ucsc.edu/cgi-bin/das/dsn/",
        "https://genome-asia.ucsc.edu/cgi-bin/das/dsn/",
        "https://genome-euro.ucsc.edu/cgi-bin/das/dsn",
    ],
    "goldenpath": [
        "https://hgdownload.cse.ucsc.edu/goldenPath",
        "https://hgdownload.soe.ucsc.edu/goldenPath",
    ],
}

available_servers = {key: input_validator.get_first_available_url(key, urls) for key, urls in server_lists.items()}


@pytest.mark.slow
def test_available_genome_pass():
    genome = "mm10"
    url_das = available_servers["das"]
    assert input_validator.is_genome_id_available_in_ucsc(genome, url_das) is True


@pytest.mark.slow
def test_available_genome_fail():
    genome = "mm12345"
    url_das = available_servers["das"]
    assert input_validator.is_genome_id_available_in_ucsc(genome, url_das) is False
