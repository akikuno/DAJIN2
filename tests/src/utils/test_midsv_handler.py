from __future__ import annotations

import pytest

from DAJIN2.utils import midsv_handler

###########################################################
# find_n_boundaries
###########################################################


@pytest.mark.parametrize(
    "midsv_tags, expected",
    [
        (["=N", "=N", "A", "B", "=N", "=N"], (1, 4)),
        (["=N", "=N", "A", "B", "A", "B"], (1, 6)),
        (["A", "B", "A", "B", "=N", "=N"], (-1, 4)),
        (["A", "B", "A", "B", "A", "B"], (-1, 6)),
    ],
)
def test_find_n_boundaries(midsv_tags, expected):
    assert midsv_handler.find_n_boundaries(midsv_tags) == expected


###########################################################
# convert midsv_tags to DNA sequence
###########################################################


@pytest.mark.parametrize(
    "sequence, expected",
    [
        ("AGGCGAacgAaccN", "AGGCGAcgtAggtN"),  # Lowercase inversion
        ("AAAAaaaaAAAA", "AAAAttttAAAA"),  # All lowercase in the middle
        ("aaaa", "tttt"),  # Only lowercase
        ("AGGCGATACC", "AGGCGATACC"),  # No lowercase
        ("a", "t"),  # Single lowercase
        ("", ""),  # Empty string
    ],
)
def test_revcomp_inversion(sequence, expected):
    assert midsv_handler._revcomp_inversion(sequence) == expected


@pytest.mark.parametrize(
    "midsv_tags, expected",
    [
        (
            ["=A", "=n", "+a|+t|+g|=a", "=c", "=g", "=A", "=a", "=c", "=c", "=N"],
            "AcgtcatnAggtN",
        ),  # Complex case with insertion and inversion
        (["=A", "=T", "=C", "=G"], "ATCG"),  # Simple match only
        (["+a|+t|+g|=a"], "tcat"),  # Single insertion with match
        (["=a", "=c", "=t"], "agt"),  # Simple inversion
        ([], ""),  # Empty midsv_tags
    ],
)
def test_convert_midsvs_to_sequence(midsv_tags, expected):
    assert midsv_handler.convert_midsvs_to_sequence(midsv_tags) == expected


###########################################################
# convert midsv_tags to cstag
###########################################################


@pytest.mark.parametrize(
    "midsv_tags, expected",
    [
        ([], []),
        (["=A", "=C", "=G"], ["=A", "=C", "=G"]),
        (["=N", "=C", "=N"], ["=N", "=C", "=N"]),
        (["=A", "+G|+G|+G|=N", "=T"], ["=A", "+G|+G|+G|=N", "=T"]),
    ],
)
def test_add_match_operator_to_n(midsv_tags, expected):
    assert midsv_handler._add_match_operator_to_n(midsv_tags) == expected


@pytest.mark.parametrize(
    "input_midsv_tags, expected",
    [
        ([], []),
        (["=ACGT"], ["=ACGT"]),
        (["=ACGT", "*GC", "=ACGT"], ["=ACGT", "*gc", "=ACGT"]),
    ],
)
def test_standardize_case(input_midsv_tags, expected):
    assert midsv_handler._standardize_case(input_midsv_tags) == expected


@pytest.mark.parametrize(
    "midsv_tags, expected",
    [
        ([], ""),
        (["=A", "=C", "=G", "=T"], "=ACGT"),
        (["=A", "*GT", "=A"], "=A*gt=A"),
        (["-C", "-G"], "-cg"),
        (["*GC", "*TA"], "*gc*ta"),
        (["+A|+A|=C", "=G"], "+aa=CG"),
        (["+A|+C|=T", "+A|+C|=G"], "+ac=T+ac=G"),
        (["=A", "=C", "+A|+C|=G", "-C", "-G", "*CG"], "=AC+ac=G-cg*cg"),
    ],
)
def test_convert_midsvs_to_cstag(midsv_tags, expected):
    assert midsv_handler.convert_midsvs_to_cstag(midsv_tags) == expected


###########################################################
# call sequence
###########################################################


@pytest.mark.parametrize(
    "cons_percentage, expected_sequence",
    [
        ([{"=A": 1.0}, {"=T": 0.9, "-T": 0.1}], "AT"),  # match
        ([{"=A": 1.0}, {"-A": 0.9, "=A": 0.1}, {"=T": 1.0}], "AT"),  # deletion
        ([{"=A": 1.0}, {"*AC": 0.9, "-A": 0.1}, {"=T": 1.0}], "ACT"),  # substitution
        ([{"=A": 1.0}, {"=a": 0.9, "-a": 0.1}, {"=T": 1.0}], "AAT"),  # inversion (not reflected in sequence...)
        ([{"=A": 1.0}, {"+G|+G|+G|=A": 0.9, "-A": 0.1}, {"=T": 1.0}], "AGGGAT"),  # insertion match
        ([{"=A": 1.0}, {"+G|+G|+G|-A": 0.9, "-A": 0.1}, {"=T": 1.0}], "AGGGT"),  # insertion deletion
        ([{"=A": 1.0}, {"+G|+G|+G|*AT": 0.9, "-A": 0.1}, {"=T": 1.0}], "AGGGTT"),  # insertion substitution
        ([{"=A": 1.0}, {"+G|+G|+G|=N": 0.9, "-A": 0.1}, {"=T": 1.0}], "AGGGNT"),  # insertion N
        ([{"=A": 1.0}], "A"),
    ],
)
def test_call_sequence(cons_percentage, expected_sequence):
    assert midsv_handler.call_sequence(cons_percentage) == expected_sequence
