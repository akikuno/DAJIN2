from __future__ import annotations

import pytest
from DAJIN2.utils.cssplits_handler import (
    find_n_boundaries,
    add_match_operator_to_n,
    convert_cssplits_to_cstag,
    standardize_case,
    call_sequence,
    get_index_of_large_deletions,
    adjust_cs_insertion,
)


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
    assert find_n_boundaries(cssplits) == expected


###########################################################
# convert cssplits to cstag
###########################################################


@pytest.mark.parametrize(
    "cssplits, expected", [([], []), (["=A", "=C", "=G"], ["=A", "=C", "=G"]), (["N", "=C", "N"], ["=N", "=C", "=N"])]
)
def test_add_match_operator_to_n(cssplits, expected):
    assert add_match_operator_to_n(cssplits) == expected


@pytest.mark.parametrize(
    "input_str, expected",
    [
        ("", ""),
        ("*AG-DEF+GHI~jkl=XYZ", "*ag-def+ghi~jkl=XYZ"),
        ("=ABCDEF", "=ABCDEF"),
        ("*abcDEF+ghi-JKL=XYZ123*AC", "*abcdef+ghi-jkl=XYZ123*ac"),
        ("*abc", "*abc"),
        ("-DEF", "-def"),
        ("+GHI", "+ghi"),
        ("=XYZ", "=XYZ"),
    ],
)
def test_standardize_case(input_str, expected):
    assert standardize_case(input_str) == expected


@pytest.mark.parametrize(
    "cssplits, expected",
    [
        ([], ""),
        (["=A", "C", "G", "T"], "=ACGT"),
        (["=A", "*GT", "=A"], "=A*gt=A"),
        (["-C", "-G"], "-cg"),
        (["*GC", "*TA"], "*gc*ta"),
        (["+A|+A|=C", "=G"], "+aa=CG"),
        (["+A|*GC|=T", "+A|+C|=G"], "+ac=T+ac=G"),
        (["+A|*GC|=T", "+A|+C|=G", "+A|+C|=G"], "+ac=T+ac=G+ac=G"),
        (["=A", "=C", "+A|+C|=G", "-C", "-G", "*CG"], "=AC+ac=G-cg*cg"),
    ],
)
def test_convert_cssplits_to_cstag(cssplits, expected):
    assert convert_cssplits_to_cstag(cssplits) == expected


###########################################################
# call sequence
###########################################################


@pytest.mark.parametrize(
    "cons_percentage_by_key, expected_sequence",
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
def test_call_sequence(cons_percentage_by_key, expected_sequence):
    assert call_sequence(cons_percentage_by_key) == expected_sequence


###########################################################
# get_index_of_large_deletions
###########################################################


@pytest.mark.parametrize(
    "cssplits, expected",
    [
        (["=T"] * 100 + ["-A"] * 300 + ["=T"] * 100 + ["-A"] * 300, set(range(100, 500))),
    ],
)
def test_get_index_of_large_deletions(cssplits, expected):
    assert get_index_of_large_deletions(cssplits) == expected


###########################################################
# reallocate_insertion_within_deletion
###########################################################


# @pytest.mark.parametrize(
#     "index, cssplits, del_range, distance, expected",
#     [
#         (0, ["-A", "=C", "=C", "=C", "=C", "=C", "=C", "=C", "=C", "=C", "=C", "-T", "=G"], 1, 10, True),
#         (11, ["-A", "=C", "=C", "=C", "=C", "=C", "=C", "=C", "=C", "=C", "=C", "-T", "=G"], 1, 10, False),
#         (0, ["-A", "=C", "=C", "=C", "=C", "=C", "=C", "=C", "=C", "=C", "=C", "=C", "-T", "=G"], 1, 10, False),
#         (0, ["-A", "=C", "=C", "=C", "=C", "=C", "*GC", "=C", "=C", "=C", "=C", "-T", "=G"], 1, 10, True),
#         (0, ["-A", "=C", "=C", "=C", "=C", "=C", "-T", "=C", "=C", "=C", "=C", "-T", "=G"], 1, 10, True),
#         (0, ["-A", "=C", "=C", "=C", "=C", "=C", "+A|=C", "=C", "=C", "=C", "=C", "-T", "=G"], 1, 10, True),
#         (0, ["-A", "-A", "=C", "-T", "-A"], 2, 10, True),
#     ],
#     ids=[
#         "-A has a 10-character match and is within the deletion cluster.",
#         "-T is outside the deletion cluster.",
#         "-A has a 11-character match and is outside the deletion cluster.",
#         "-A is determined to be reset and within the deletion cluster due to a point mutation at the 6th character.",
#         "-A is determined to be reset and within the deletion cluster due to a deletion at the 6th character.",
#         "-A is determined to be reset and within the deletion cluster due to an insertion at the 6th character.",
#         "-A is within the deletion cluster.",
#     ],
# )
# def test_is_within_deletion(index: int, cssplits: list[str], del_range: int, distance: int, expected: bool):
#     assert is_within_deletion(index, cssplits, del_range, distance) == expected


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
    assert adjust_cs_insertion(cs) == expected


# @pytest.mark.parametrize(
#     "input_str, expected_output",
#     [
#         ("-A,-A,-A,=C,=C,=C,-T,-T,-T,=G", "-A,-A,-A,-C,-C,-C,-T,-T,-T,+C|+C|+C|=G"),
#         ("-A,-A,-A,=C,=C,=C,=C,-T,-T,-T", "-A,-A,-A,=C,=C,=C,=C,-T,-T,-T"),
#         ("-A,-A,-A,N,=C,n,-T,-T,-T,=G", "-A,-A,-A,N,-C,n,-T,-T,-T,+N|+C|+n|=G"),
#         ("-A,-A,-A,=C,+T|+T|=C,=C,-T,-T,-T,=G", "-A,-A,-A,-C,-C,-C,-T,-T,-T,+C|+T|+T|+C|+C|=G"),
#         ("-A,-A,-A,=C,+T|+T|*CG,=C,-T,-T,-T,=G", "-A,-A,-A,-C,-C,-C,-T,-T,-T,+C|+T|+T|+G|+C|=G"),
#         ("-G,-G,-C,=A,=C,=C,*CA,=A,-T,-T,*AC", "-G,-G,-C,=A,=C,=C,*CA,=A,-T,-T,*AC"),
#     ],
#     ids=[
#         "insertion within deletion",
#         "4-character match",
#         "N and n",
#         "Insertion",
#         "Insertion followed by substitution",
#         "Should not be adjusted",
#     ],
# )
# def test_reallocate_insertion_within_deletion(input_str: str, expected_output: str):
#     assert reallocate_insertion_within_deletion(input_str, del_range=3, distance=3) == expected_output
