from __future__ import annotations

import pytest

from DAJIN2.utils import cssplits_handler

###########################################################
# find_n_boundaries
###########################################################


@pytest.mark.parametrize(
    "cssplits, expected",
    [
        (["N", "N", "A", "B", "N", "N"], (1, 4)),
        (["N", "N", "A", "B", "A", "B"], (1, 6)),
        (["A", "B", "A", "B", "N", "N"], (-1, 4)),
        (["A", "B", "A", "B", "A", "B"], (-1, 6)),
    ],
)
def test_find_n_boundaries(cssplits, expected):
    assert cssplits_handler.find_n_boundaries(cssplits) == expected


###########################################################
# convert cssplits to cstag
###########################################################


@pytest.mark.parametrize(
    "cssplits, expected",
    [
        ([], []),
        (["=A", "=C", "=G"], ["=A", "=C", "=G"]),
        (["N", "=C", "N"], ["=N", "=C", "=N"]),
        (["=A", "+G|+G|+G|N", "=T"], ["=A", "+G|+G|+G|=N", "=T"]),
    ],
)
def test_add_match_operator_to_n(cssplits, expected):
    assert cssplits_handler._add_match_operator_to_n(cssplits) == expected


@pytest.mark.parametrize(
    "input_cssplits, expected",
    [
        ([], []),
        (["=ACGT"], ["=ACGT"]),
        (["=ACGT", "*GC", "=ACGT"], ["=ACGT", "*gc", "=ACGT"]),
    ],
)
def test_standardize_case(input_cssplits, expected):
    assert cssplits_handler._standardize_case(input_cssplits) == expected


@pytest.mark.parametrize(
    "cssplits, expected",
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
def test_convert_cssplits_to_cstag(cssplits, expected):
    assert cssplits_handler.convert_cssplits_to_cstag(cssplits) == expected


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
        ([{"=A": 1.0}, {"+G|+G|+G|N": 0.9, "-A": 0.1}, {"=T": 1.0}], "AGGGNT"),  # insertion N
        ([{"=A": 1.0}], "A"),
    ],
)
def test_call_sequence(cons_percentage, expected_sequence):
    assert cssplits_handler.call_sequence(cons_percentage) == expected_sequence


###########################################################
# get_index_of_large_deletions
###########################################################


@pytest.mark.parametrize(
    "cssplits, expected",
    [
        (["=T"] * 100 + ["-A"] * 300 + ["=T"] * 100, set(range(100, 400))),
        (
            ["=T"] * 100 + ["-A"] * 300 + ["=T"] * 10 + ["-A"] * 300 + ["=T"] * 100,
            set(range(100, 400)) | set(range(410, 710)),
        ),
    ],
)
def test_get_index_of_large_deletions(cssplits, expected):
    assert cssplits_handler._get_index_of_large_deletions(cssplits) == expected


###########################################################
# reallocate_insertion_within_deletion
###########################################################


@pytest.mark.parametrize(
    "cs, expected",
    [
        ("=G", "+G|"),
        ("*AG", "+G|"),
        ("+A|+C|+G|=T", "+A|+C|+G|+T|"),
        ("N", "+N|"),
        ("-G", None),
    ],
)
def test_adjust_cs_insertion(cs: str, expected: str):
    assert cssplits_handler._adjust_cs_insertion(cs) == expected


@pytest.mark.parametrize(
    "cssplits, expected",
    [
        (
            ["=T"] * 100 + ["-A"] * 300 + ["*TA"] * 10 + ["-A"] * 300 + ["=T"] * 100,
            ["=T"] * 100
            + ["-A"] * 300
            + ["-T"] * 10
            + ["-A"] * 300
            + ["+A|+A|+A|+A|+A|+A|+A|+A|+A|+A|=T"]
            + ["=T"] * 99,
        ),
        (
            ["=T"] * 100 + ["-A"] * 150 + ["=T"] * 10 + ["-A"] * 150 + ["=T"] * 100,
            ["=T"] * 100 + ["-A"] * 150 + ["=T"] * 10 + ["-A"] * 150 + ["=T"] * 100,
        ),
        (
            ["=T"] * 100
            + ["-A"] * 100
            + ["*TA"] * 10
            + ["-A"] * 100
            + ["=T"] * 10
            + ["-A"] * 100
            + ["*TA"] * 10
            + ["-A"] * 100
            + ["=T"] * 100,
            ["=T"] * 100
            + ["-A"] * 100
            + ["-T"] * 10
            + ["-A"] * 100
            + ["+A|+A|+A|+A|+A|+A|+A|+A|+A|+A|=T"]
            + ["=T"] * 9
            + ["-A"] * 100
            + ["-T"] * 10
            + ["-A"] * 100
            + ["+A|+A|+A|+A|+A|+A|+A|+A|+A|+A|=T"]
            + ["=T"] * 99,
        ),
    ],
    ids=[
        "insertion within deletion",
        "matched region within deletion",
        "insertions within deletion and matched region",
    ],
)
def test_reallocate_insertion_within_deletion(cssplits: str, expected: str):
    assert cssplits_handler.reallocate_insertion_within_deletion(cssplits) == expected
