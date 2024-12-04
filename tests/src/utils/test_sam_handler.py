from __future__ import annotations

from pathlib import Path

import midsv

from DAJIN2.utils.sam_handler import (
    calculate_alignment_length,
    is_header,
    is_mapped,
    is_overlapping,
    remove_microhomology,
    remove_overlapped_reads,
    revcomp_sam,
    reverse_flag,
    split_cigar,
)


def test_split_cigar():
    assert split_cigar("10M2D3M") == ["10M", "2D", "3M"]
    assert split_cigar("5M1I5M") == ["5M", "1I", "5M"]
    assert split_cigar("3M1D1I3M") == ["3M", "1D", "1I", "3M"]
    assert split_cigar("10M") == ["10M"]
    assert split_cigar("") == []
    assert split_cigar("3S10M2D3M") == ["3S", "10M", "2D", "3M"]
    assert split_cigar("5M1I5M3H") == ["5M", "1I", "5M", "3H"]
    assert split_cigar("3S3M1D1I3M3H") == ["3S", "3M", "1D", "1I", "3M", "3H"]


def test_calculate_alignment_length():
    assert calculate_alignment_length("10M2D3M") == 15
    assert calculate_alignment_length("5M1I5M") == 10
    assert calculate_alignment_length("3M1D1I3M") == 7
    assert calculate_alignment_length("3S10M2D3M") == 15
    assert calculate_alignment_length("5M1I5M3H") == 10
    assert calculate_alignment_length("3S3M1D1I3M3H") == 7
    assert calculate_alignment_length("10M") == 10
    assert calculate_alignment_length("") == 0


###########################################################
# remove_overlapped_reads
###########################################################


def test_is_header():
    assert is_header(["@HD"])
    assert not is_header(["read1", "flag", "rname", "3", "mapq", "4M", "rnext", "pnext", "tlen", "ACTG"])


def test_is_mapped():
    assert not is_mapped(["@HD"])
    assert is_mapped(["read1", "flag", "rname", "3", "mapq", "4M", "rnext", "pnext", "tlen", "ACTG"])


def test_is_overlapping():
    assert is_overlapping(
        ["read1", "flag", "rname", "3", "mapq", "4M", "rnext", "pnext", "tlen", "ACTG"],
        ["read1", "flag", "rname", "5", "mapq", "4M", "rnext", "pnext", "tlen", "ACTG"],
    )
    assert not is_overlapping(
        ["read1", "flag", "rname", "3", "mapq", "4M", "rnext", "pnext", "tlen", "ACTG"],
        ["read1", "flag", "rname", "8", "mapq", "4M", "rnext", "pnext", "tlen", "ACTG"],
    )


def test_remove_overlapped_reads_simulation():
    # Case with only headers
    sam = [["@header1"], ["@header2"]]
    assert remove_overlapped_reads(sam) == sam

    # Case with headers and content, but the content is unmapped ("*")
    sam = [["@header1"], ["read1", "flag", "rname", "3", "mapq", "*", "rnext", "pnext", "tlen", "*"]]
    assert remove_overlapped_reads(sam) == [["@header1"]]

    # Case with the same read name mapped to different positions
    sam = [
        ["@header1"],
        ["read1", "flag", "rname", "3", "mapq", "4M", "rnext", "pnext", "tlen", "ACTG"],
        ["read1", "flag", "rname", "10", "mapq", "4M", "rnext", "pnext", "tlen", "ACTG"],
    ]
    assert remove_overlapped_reads(sam) == sam

    # Case with the same read name mapped to overlapping positions
    sam = [
        ["@header1"],
        ["read1", "flag", "rname", "3", "mapq", "4M", "rnext", "pnext", "tlen", "ACTG"],
        ["read1", "flag", "rname", "4", "mapq", "4M", "rnext", "pnext", "tlen", "ACTG"],
    ]
    assert remove_overlapped_reads(sam) == [
        ["@header1"],
    ]


def test_remove_overlapped_reads_real_sample():
    """
    # Overlapped reads in tyr albino:
        - 01f0e046-b6dc-46f3-ad5b-96ed91484778
    # Non-overlapped reads in tyr albino:
        - a224e9ca-d634-4490-bf77-abd9d5cbd7bc
    """
    sam = midsv.read_sam("tests/data/report/report_bam/remove_overlap.sam")
    sam = list(sam)
    test = remove_overlapped_reads(sam)
    answer = midsv.read_sam("tests/data/report/report_bam/answer.sam")
    answer = list(answer)
    assert test == answer


###########################################################
# remove_microholomogy
###########################################################


def test_remove_microhomology_ACGT():
    sam = midsv.read_sam("tests/data/report/microhomology/microhomology-ACGT.sam")
    sam = list(sam)
    test = remove_microhomology(sam)
    answer = Path("tests/data/report/microhomology/microhomology-ACGT-answer.txt").read_text()
    answer = eval(answer)
    assert test == answer


def test_remove_microhomology_insertion():
    sam = midsv.read_sam("tests/data/report/microhomology/microhomology-insertion.sam")
    sam = list(sam)
    test = remove_microhomology(sam)
    answer = Path("tests/data/report/microhomology/microhomology-insertion-answer.txt").read_text()
    answer = eval(answer)
    assert test == answer


def test_remove_microhomology_deletion():
    sam = midsv.read_sam("tests/data/report/microhomology/microhomology-deletion.sam")
    sam = list(sam)
    test = remove_microhomology(sam)
    answer = Path("tests/data/report/microhomology/microhomology-deletion-answer.txt").read_text()
    answer = eval(answer)
    assert test == answer


def test_remove_microhomology_overlapped_softclip():
    sam = midsv.read_sam("tests/data/report/report_bam/cables2_barcode34_microhomology.sam")
    sam = list(sam)
    test = remove_microhomology(sam)
    answer = list(midsv.read_sam("tests/data/report/report_bam/answer_cables2_barcode34_microhomology.sam"))
    assert test == answer


def test_remove_microhomology_real_singe_read():
    seq_id = "0ef147016ef0"
    sam = midsv.read_sam("tests/data/report/report_bam/barcode54_allele2_before.sam")
    sam = list(sam)
    sam_headers = [s for s in sam if s[0].startswith("@")]
    sam_contents = [s for s in sam if seq_id in s[0] and s[9] != "*"]
    sam = sam_headers + sam_contents
    test = remove_microhomology(sam)
    answer = midsv.read_sam(f"tests/data/report/report_bam/barcode54_allele2_{seq_id}_after.sam")
    answer = list(answer)
    assert test == answer


def test_remove_microhomology_real_all_500_reads():
    sam = midsv.read_sam("tests/data/report/report_bam/barcode54_allele2_before.sam")
    sam = list(sam)
    test = remove_microhomology(sam)
    answer = midsv.read_sam("tests/data/report/report_bam/barcode54_allele2_after.sam")
    answer = list(answer)
    assert test == answer


###########################################################
# revcomp_sam
###########################################################


def test_reverse_flag():
    # Forward strand (0) to reverse strand (16)
    assert reverse_flag(0) == 16

    # Reverse strand (16) to forward strand (0)
    assert reverse_flag(16) == 0

    # Forward strand with supplementary alignment (2048) to reverse strand with supplementary alignment (2064)
    assert reverse_flag(2048) == 2064

    # Reverse strand with supplementary alignment (2064) to forward strand with supplementary alignment (2048)
    assert reverse_flag(2064) == 2048

    # Other flags should remain unchanged when flipped
    # For example: 73 (01001001) -> 89 (01011001) -> 73 (01001001)
    assert reverse_flag(reverse_flag(73)) == 73


def test_revcomp_sam():
    # Single entry in the SAM file with different quality scores
    sam_contents = [["read1", "0", "chr1", "1", "60", "4M", "*", "0", "0", "ATGC", "!@#$"]]
    genome_end = 5
    assert revcomp_sam(sam_contents, genome_end) == [
        ["read1", "16", "chr1", "2", "60", "4M", "*", "0", "0", "GCAT", "$#@!"]
    ]

    # Multiple entries in the SAM file with different flags and quality scores
    sam_contents = [
        ["read1", "0", "chr1", "1", "60", "2M2D2M", "*", "0", "0", "ATGC", "!@#$"],
        ["read2", "16", "chr1", "1", "60", "2M2I2M", "*", "0", "0", "GCTA", "abcd"],
    ]
    genome_end = 6
    assert revcomp_sam(sam_contents, genome_end) == [
        ["read1", "16", "chr1", "1", "60", "2M2D2M", "*", "0", "0", "GCAT", "$#@!"],
        ["read2", "0", "chr1", "3", "60", "2M2I2M", "*", "0", "0", "TAGC", "dcba"],
    ]

    # Flag not in flag_map, with different quality scores
    sam_contents = [["read1", "4", "chr1", "1", "60", "4M", "*", "0", "0", "ATGC", "1234"]]
    genome_end = 5
    assert revcomp_sam(sam_contents, genome_end) == [
        ["read1", "20", "chr1", "2", "60", "4M", "*", "0", "0", "GCAT", "4321"]
    ]
