import pytest
from DAJIN2.utils.cssplits_handler import (
    find_n_boundaries,
    add_match_operator_to_n,
    concatenate_cssplits,
    standardize_case,
    call_sequence,
    is_start_of_deletion,
    is_within_deletion,
    adjust_cs_insertion,
    detect_insertion_within_deletion,
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


def test_add_match_operator_to_n_empty_input():
    assert add_match_operator_to_n([]) == []


def test_add_match_operator_to_n_no_n_starting_strings():
    cssplits = ["=A", "=C", "=G"]
    assert add_match_operator_to_n(cssplits) == cssplits


def test_add_match_operator_to_n_with_n_starting_strings():
    cssplits = ["N", "=C", "N"]
    expected_output = ["=N", "=C", "=N"]
    assert add_match_operator_to_n(cssplits) == expected_output


def test_concatenate_cssplits_empty_input():
    assert concatenate_cssplits([]) == ""


def test_concatenate_cssplits_single_element():
    assert concatenate_cssplits(["=ACGT"]) == "=ACGT"


def test_concatenate_cssplits_same_symbols():
    cssplits = ["=AC", "=GT", "-CC", "-GG"]
    assert concatenate_cssplits(cssplits) == "=ACGT-CCGG"


def test_concatenate_cssplits_plus_symbol():
    cssplits = ["+A|*GC|=T", "+A|+C|=G"]
    assert concatenate_cssplits(cssplits) == "+AC=T+AC=G"


def test_concatenate_cssplits_plus_symbol_torio():
    cssplits = ["+A|*GC|=T", "+A|+C|=G", "+A|+C|=G"]
    assert concatenate_cssplits(cssplits) == "+AC=T+AC=G+AC=G"


def test_concatenate_cssplits_mixed():
    cssplits = ["=AC", "=GT", "+A|+C|=G", "-CC", "-GG"]
    assert concatenate_cssplits(cssplits) == "=ACGT+AC=G-CCGG"


def test_standardize_case_empty_input():
    assert standardize_case("") == ""


def test_standardize_case_special_characters():
    assert standardize_case("*AG-DEF+GHI~jkl=XYZ") == "*ag-def+ghi~jkl=XYZ"


def test_standardize_case_no_special_characters():
    assert standardize_case("=ABCDEF") == "=ABCDEF"


def test_standardize_case_mixed():
    assert standardize_case("*abcDEF+ghi-JKL=XYZ123*AC") == "*abcdef+ghi-jkl=XYZ123*ac"


@pytest.mark.parametrize(
    "input, expected",
    [
        ("*abc", "*abc"),
        ("-DEF", "-def"),
        ("+GHI", "+ghi"),
        ("=XYZ", "=XYZ"),
    ],
)
def test_standardize_case_parametrized(input, expected):
    assert standardize_case(input) == expected


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
# detect_insertion_within_deletion
###########################################################


@pytest.mark.parametrize(
    "i, cssplits, start_del_count, expected",
    [
        (0, ["-A", "-A", "-A"], 3, True),
        (0, ["-A", "-A"], 3, False),
        (1, ["-A", "-A", "-A"], 3, False),
    ],
)
def test_is_start_of_deletion(i: int, cssplits: list[str], start_del_count: int, expected: bool):
    assert is_start_of_deletion(i, cssplits, start_del_count) == expected


@pytest.mark.parametrize(
    "index, cssplits, expected",
    [
        (0, ["-A", "=C", "=C", "=C", "=C", "=C", "=C", "=C", "=C", "=C", "=C", "-T", "=G"], True),
        (11, ["-A", "=C", "=C", "=C", "=C", "=C", "=C", "=C", "=C", "=C", "=C", "-T", "=G"], False),
        (0, ["-A", "=C", "=C", "=C", "=C", "=C", "=C", "=C", "=C", "=C", "=C", "=C", "-T", "=G"], False),
        (0, ["-A", "=C", "=C", "=C", "=C", "=C", "*GC", "=C", "=C", "=C", "=C", "-T", "=G"], True),
        (0, ["-A", "=C", "=C", "=C", "=C", "=C", "-T", "=C", "=C", "=C", "=C", "-T", "=G"], True),
        (0, ["-A", "=C", "=C", "=C", "=C", "=C", "+A|=C", "=C", "=C", "=C", "=C", "-T", "=G"], True),
        (0, ["-A", "=C", "=C", "-T", "=G"], True),
    ],
    ids=[
        "-A has a 10-character match and is within the deletion cluster.",
        "-T is outside the deletion cluster.",
        "-A has a 11-character match and is outside the deletion cluster.",
        "-A is determined to be reset and within the deletion cluster due to a point mutation at the 6th character.",
        "-A is determined to be reset and within the deletion cluster due to a deletion at the 6th character.",
        "-A is determined to be reset and within the deletion cluster due to an insertion at the 6th character.",
        "-A is within the deletion cluster.",
    ],
)
def test_is_within_deletion(index: int, cssplits: list[str], expected: bool):
    assert is_within_deletion(index, cssplits) == expected


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


@pytest.mark.parametrize(
    "input_str, expected_output",
    [
        ("=A,-A,-A,-A,=C,=C,=C,-T,=G", "=A,-A,-A,-A,-C,-C,-C,-T,+C|+C|+C|=G"),
        (
            "-A,-A,-A,=C,=C,=C,=C,=C,=C,=C,=C,=C,=C,-T,=G",
            "-A,-A,-A,-C,-C,-C,-C,-C,-C,-C,-C,-C,-C,-T,+C|+C|+C|+C|+C|+C|+C|+C|+C|+C|=G",
        ),
        ("-A,-A,-A,=C,=C,=C,=C,=C,=C,=C,=C,=C,=C,=C,-T,=G", "-A,-A,-A,=C,=C,=C,=C,=C,=C,=C,=C,=C,=C,=C,-T,=G"),
        ("=A,-A,-A,-A,N,=C,n,-T,=G", "=A,-A,-A,-A,-N,-C,-n,-T,+N|+C|+n|=G"),
        ("=A,-A,-A,-A,=C,+T|+T|=C,=C,-T,=G", "=A,-A,-A,-A,-C,-C,-C,-T,+C|+T|+T|+C|+C|=G"),
        ("=A,-A,-A,-A,=C,+T|+T|*CG,=C,-T,=G", "=A,-A,-A,-A,-C,-C,-C,-T,+C|+T|+T|+G|+C|=G"),
    ],
    ids=[
        "insertion within deletion",
        "10-character match",
        "11-character match",
        "N and n",
        "Insertion",
        "Insertion followed by substitution",
    ],
)
def test_detect_insertion_within_deletion(input_str: str, expected_output: str):
    assert detect_insertion_within_deletion(input_str) == expected_output
