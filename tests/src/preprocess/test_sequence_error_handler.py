from DAJIN2.core.preprocess.error_correction.sequence_error_fastq import (
    extract_qname_from_fastq_identifier,
    format_fastq_record,
    split_fastq_records_by_qname,
)


def test_extract_qname_from_fastq_identifier_ignores_metadata():
    assert extract_qname_from_fastq_identifier("@read1 runid=abc") == "read1"
    assert extract_qname_from_fastq_identifier("read2") == "read2"


def test_split_fastq_records_by_qname_separates_passed_and_error_reads():
    fastq_records = [
        {"identifier": "@read1 runid=abc", "sequence": "AAAA", "separator": "+", "quality": "!!!!"},
        {"identifier": "@read2", "sequence": "CCCC", "separator": "+", "quality": "####"},
        {"identifier": "@read3", "sequence": "GGGG", "separator": "+", "quality": "$$$$"},
    ]

    passed, error = split_fastq_records_by_qname(fastq_records, {"read1", "read3"})

    assert [record["identifier"] for record in passed] == ["@read1 runid=abc", "@read3"]
    assert [record["identifier"] for record in error] == ["@read2"]


def test_format_fastq_record_writes_four_line_record():
    read = {"identifier": "@read1", "sequence": "ACGT", "separator": "+", "quality": "!!!!"}

    assert format_fastq_record(read) == "@read1\nACGT\n+\n!!!!\n"
