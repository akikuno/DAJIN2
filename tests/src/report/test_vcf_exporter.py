from __future__ import annotations

import sys
import types
from pathlib import Path

for module_name in ("pysam", "cstag", "mappy", "wslPath"):
    try:
        __import__(module_name)
    except ModuleNotFoundError:
        sys.modules[module_name] = types.SimpleNamespace()

from src.DAJIN2.core import report
from src.DAJIN2.utils.midsv_handler import revcomp_midsvs


def test_midsv_to_vcf_records_insertion_uses_left_anchor_position():
    midsv_tags = ["=A", "+G|+T|=A", "=G"]
    records = report.vcf_exporter._midsv_to_vcf_records(midsv_tags, "chr1", 0, "allele01", 50)

    assert len(records) == 1
    record = records[0]
    assert record["POS"] == 1
    assert record["REF"] == "A"
    assert record["ALT"] == "AGT"
    assert record["INFO"]["SVLEN"] == 2


def test_midsv_to_vcf_records_deletion_uses_deleted_base_position_and_end():
    midsv_tags = ["=A", "-G", "-T", "=G"]
    records = report.vcf_exporter._midsv_to_vcf_records(midsv_tags, "chr1", 0, "allele01", 50)

    assert len(records) == 1
    record = records[0]
    assert record["POS"] == 2
    assert record["REF"] == "G"
    assert record["ALT"] == "<DEL>"
    assert record["INFO"]["SVLEN"] == -2
    assert record["INFO"]["SEQ"] == "GT"
    assert record["INFO"]["END"] == 3


def test_midsv_to_vcf_records_deletion_with_offset_sets_expected_coordinates():
    midsv_tags = ["=A", "-G", "=T"]
    records = report.vcf_exporter._midsv_to_vcf_records(midsv_tags, "chr7", 100, "allele01", 50)

    assert len(records) == 1
    record = records[0]
    assert record["POS"] == 102
    assert record["REF"] == "G"
    assert record["INFO"]["END"] == 102


def test_midsv_to_vcf_records_insertion_after_revcomp_keeps_left_anchor_logic():
    midsv_tags = ["=A", "+G|+T|=A", "=G"]
    midsv_tags_revcomp = revcomp_midsvs(midsv_tags)
    records = report.vcf_exporter._midsv_to_vcf_records(midsv_tags_revcomp, "chr1", 0, "allele01", 50)

    assert len(records) == 1
    record = records[0]
    assert record["POS"] == 2
    assert record["REF"] == "T"
    assert record["ALT"] == "TAC"


def test_export_to_vcf_skips_false_positive_when_consensus_matches_control(tmp_path, monkeypatch):
    monkeypatch.setattr(report.vcf_exporter, "_load_control_sequence", lambda _tempdir, _sample_name: "ACGT")
    monkeypatch.setattr(report.vcf_exporter, "convert_midsvs_to_sequence", lambda _tags: "ACGT")
    monkeypatch.setattr(
        report.vcf_exporter, "_align_sequence_to_control", lambda _control, _query: ["=A", "=C", "=G", "=T"]
    )

    report.vcf_exporter.export_to_vcf(
        tmp_path,
        "sample",
        {"chrom": "chr1", "start": 100, "strand": "+"},
        {"allele02_control_intact_35.0%": ["=A"]},
    )

    path_vcf = Path(tmp_path, "report", "VCF", "sample", "sample_allele02_control_intact_35.0%.vcf")
    contents = path_vcf.read_text(encoding="utf-8").strip().splitlines()
    assert contents == ["##fileformat=VCFv4.3", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"]


def test_export_to_vcf_emits_multiple_insertions_from_consensus_alignment(tmp_path, monkeypatch):
    monkeypatch.setattr(report.vcf_exporter, "_load_control_sequence", lambda _tempdir, _sample_name: "ACGT")
    monkeypatch.setattr(report.vcf_exporter, "convert_midsvs_to_sequence", lambda _tags: "ATCGGT")
    monkeypatch.setattr(
        report.vcf_exporter,
        "_align_sequence_to_control",
        lambda _control, _query: ["=A", "+T|=A", "=C", "+G|+G|=C", "=G", "=T"],
    )

    report.vcf_exporter.export_to_vcf(
        tmp_path,
        "sample",
        {"chrom": "chr1", "start": 100, "strand": "+"},
        {"allele01_flox_intact_53.021%": ["=A"]},
    )

    path_vcf = Path(tmp_path, "report", "VCF", "sample", "sample_allele01_flox_intact_53.021%.vcf")
    records = path_vcf.read_text(encoding="utf-8").strip().splitlines()[2:]

    assert len(records) == 2
    assert records[0].startswith("chr1\t100\tallele01_flox_intact_53.021%25\tA\tAT")
    assert "SVLEN=1" in records[0]
    assert records[1].startswith("chr1\t102\tallele01_flox_intact_53.021%25\tC\tCGG")
    assert "SVLEN=2" in records[1]


def test_convert_records_to_genomic_coordinates_minus_strand_sanitizes_inserted_sequence():
    records = [
        {
            "CHROM": "control",
            "POS": 10,
            "ID": "allele01",
            "REF": "A",
            "ALT": "AAC",
            "INFO": {"TYPE": "INS", "SVLEN": 2, "SEQ": "=AC", "QNAME": "allele01"},
        }
    ]

    converted = report.vcf_exporter._convert_records_to_genomic_coordinates(
        records,
        {"chrom": "chr2", "start": 100, "end": 200, "strand": "-"},
    )

    assert len(converted) == 1
    record = converted[0]
    assert record["POS"] == 191
    assert record["REF"] == "T"
    assert record["ALT"] == "TGT"
    assert record["INFO"]["SEQ"] == "GT"
