import pytest
from DAJIN2.utils.cssplits_handler import (
    find_n_boundaries,
    add_match_operator_to_n,
    concatenate_cssplits,
    standardize_case,
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
