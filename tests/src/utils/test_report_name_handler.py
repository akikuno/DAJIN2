from __future__ import annotations

from DAJIN2.utils.report_name_handler import build_report_filename


def test_build_report_filename_compacts_long_consensus_label():
    long_allele_name = "deletion_" + "A" * 400
    label = f"allele01|{long_allele_name}|indels|33.884%"

    got = build_report_filename(label, ".fasta", sample_name="sample")

    assert got == "sample_allele01_33.884%.fasta"


def test_build_report_filename_keeps_short_name():
    label = "allele01|control|intact|66.116%"
    got = build_report_filename(label, ".html", sample_name="sample")
    assert got == "sample_allele01_control_intact_66.116%.html"
