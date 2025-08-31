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
        expected = (
            f"{test} requires extensions either .fastq, .fastq.gz, .fq, .fq.gz, .fasta, .fasta.gz, .fa, .fa.gz, or .bam"
        )
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
    "gggenome": [
        "https://gggenome.dbcls.jp/",
    ],
    "goldenpath": [
        "https://hgdownload.cse.ucsc.edu/goldenPath",
        "https://hgdownload.soe.ucsc.edu/goldenPath",
    ],
}

available_servers = {key: input_validator.get_first_available_url(key, urls) for key, urls in server_lists.items()}

###############################################################################
# validate BED file
###############################################################################


@pytest.fixture
def noop_validate_exists(monkeypatch):
    """Make validate_file_existence a no-op by default."""
    monkeypatch.setattr(input_validator, "validate_file_existence", lambda p: None)


@pytest.mark.parametrize(
    "fname",
    [
        "ok.bed",
        "ok.bed.gz",
        "OK.BED",
        "OK.BED.GZ",
        "sample.region.bed",
        "sample.region.bed.gz",
    ],
)
def test_validate_bed_file_valid_extensions(tmp_path, monkeypatch, noop_validate_exists, fname):
    # The function only checks extension and calls validate_file_existence;
    # it does not read file content. Path existence is simulated by the no-op.
    bed_path = tmp_path / fname
    # Even though validate_file_existence is a no-op, create the file to be realistic
    bed_path.write_text("# dummy\n", encoding="utf-8")

    assert input_validator.validate_bed_file(str(bed_path)) is None


@pytest.mark.parametrize(
    "fname",
    [
        "bad.txt",
        "bad.bed.bgz",
        "bad.gz",
        "badbed",  # no suffix
        "bad.bed.zip",
        "bad.bedgz",  # looks similar but not .bed.gz
        "bad.bed.GZ.bak",
    ],
)
def test_validate_bed_file_invalid_extensions(tmp_path, monkeypatch, noop_validate_exists, fname):
    bed_path = tmp_path / fname
    bed_path.write_text("# dummy\n", encoding="utf-8")

    with pytest.raises(ValueError) as ei:
        input_validator.validate_bed_file(str(bed_path))
    msg = str(ei.value)
    assert "BED file must have .bed or .bed.gz extension" in msg
    assert str(bed_path) in msg


def test_validate_bed_file_wraps_bederror(monkeypatch, tmp_path):
    # Simulate a format error raised by validate_file_existence as BEDError
    err = input_validator.BEDError("malformed header")
    monkeypatch.setattr(input_validator, "validate_file_existence", lambda p: (_ for _ in ()).throw(err))

    bed_path = tmp_path / "looks_ok.bed"
    bed_path.write_text("# dummy\n", encoding="utf-8")

    with pytest.raises(ValueError) as ei:
        input_validator.validate_bed_file(str(bed_path))
    assert "Invalid BED file format: malformed header" in str(ei.value)


@pytest.mark.parametrize("exc_type", [FileNotFoundError, PermissionError, OSError, RuntimeError])
def test_validate_bed_file_wraps_other_exceptions(monkeypatch, tmp_path, exc_type):
    # Any non-BEDError exception should be wrapped with a generic message
    def raise_exc(_):
        raise exc_type("boom")

    monkeypatch.setattr(input_validator, "validate_file_existence", raise_exc)

    bed_path = tmp_path / "missing.bed"
    # Do not create the file; existence is checked by the patched function

    with pytest.raises(ValueError) as ei:
        input_validator.validate_bed_file(str(bed_path))
    msg = str(ei.value)
    assert f"Error processing BED file {bed_path}" in msg
    assert "boom" in msg


def test_validate_bed_file_accepts_double_suffix_recognition(tmp_path, monkeypatch, noop_validate_exists):
    # Ensure logic uses all suffixes (e.g., ".bed.gz") rather than only .suffix
    bed_path = tmp_path / "double.suffix.name.bed.gz"
    bed_path.write_text("# dummy\n", encoding="utf-8")
    assert input_validator.validate_bed_file(str(bed_path)) is None
