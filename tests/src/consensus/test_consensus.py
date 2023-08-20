import pytest

from src.DAJIN2.core.consensus.consensus import (
    _remove_nonconsecutive_n,
    _update_percentage,
    _call_percentage,
    _process_base,
    _call_sequence,
)


@pytest.mark.parametrize(
    "cons_percentage, expected_output",
    [
        # not remove
        (
            [{"N": 100}, {"N": 100}, {"A": 50, "C": 50}, {"N": 100}],
            [{"N": 100}, {"N": 100}, {"A": 50, "C": 50}, {"N": 100}],
        ),
        ([{"N": 80, "A": 20}, {"A": 50, "C": 50}], [{"N": 80, "A": 20}, {"A": 50, "C": 50}]),
        ([{"A": 50, "C": 50}, {"N": 80, "A": 20}], [{"A": 50, "C": 50}, {"N": 80, "A": 20}]),
        ([{"A": 50, "C": 50}], [{"A": 50, "C": 50}]),
        ([{"N": 100}], [{"N": 100}]),
        ([], []),
        # remove non consecutive N
        ([{"N": 100}, {"N": 100}, {"A": 60, "N": 40}, {"N": 100}], [{"N": 100}, {"N": 100}, {"A": 60}, {"N": 100}]),
        ([{"N": 20, "A": 80}, {"A": 50, "C": 50}], [{"A": 80}, {"A": 50, "C": 50}]),
        ([{"A": 50, "C": 50}, {"N": 20, "A": 80}], [{"A": 50, "C": 50}, {"A": 80}]),
        ([{"A": 80, "N": 20}], [{"A": 80}]),
    ],
)
def test_remove_nonconsecutive_n(cons_percentage, expected_output):
    result = _remove_nonconsecutive_n(cons_percentage)
    assert result == expected_output


def test_update_percentage():
    cons_percentage = [{"A": 25, "C": 25, "SEQERROR": 50}, {"SEQERROR": 100}]
    expected_output = [{"A": 50, "C": 50}, {"N": 100}]
    assert _update_percentage(cons_percentage) == expected_output


def test_call_percentage():
    cssplits = [["+A", "-T", "C", "A", "T"], ["+A", "=T", "C", "*AT", "*AT"]]
    mutation_loci = [{"+", "-"}, {"-"}, {}, {}, {"*"}]
    expected_output = [{"+A": 100.0}, {"-T": 50.0, "=T": 50.0}, {"C": 100.0}, {"A": 100.0}, {"T": 50.0, "*AT": 50.0}]
    assert _call_percentage(cssplits, mutation_loci) == expected_output


###########################################################
# call sequence
###########################################################


# Example test cases for _process_base function
@pytest.mark.parametrize(
    "cons, expected_output",
    [
        ("=A", "A"),
        ("-A", ""),
        ("*AC", "C"),
        ("+G|+G|+G|=A", "GGGA"),
        ("+G|+G|+G|-A", "GGG"),
        ("+G|+G|+G|*AT", "GGGT"),
        ("+G|+G|+G|N", "GGGN"),
        ("A", "A"),
        ("N", "N"),
        ("", ""),
    ],
)
def test_process_base(cons, expected_output):
    result = _process_base(cons)
    assert result == expected_output


@pytest.mark.parametrize(
    "cons_percentage_by_key, expected_sequence",
    [
        ([{"=A": 1.0}, {"=T": 0.9, "-T": 0.1}], "AT"),  # match
        ([{"=A": 1.0}, {"-A": 0.9, "=A": 0.1}, {"=T": 1.0}], "AT"),  # deletion
        ([{"=A": 1.0}, {"*AC": 0.9, "-A": 0.1}, {"=T": 1.0}], "ACT"),  # substitution
        ([{"=A": 1.0}, {"=a": 0.9, "-a": 0.1}, {"=T": 1.0}], "AaT"),  # inversion
        ([{"=A": 1.0}, {"+G|+G|+G|=A": 0.9, "-A": 0.1}, {"=T": 1.0}], "AGGGAT"),  # insertion match
        ([{"=A": 1.0}, {"+G|+G|+G|-A": 0.9, "-A": 0.1}, {"=T": 1.0}], "AGGGT"),  # insertion deletion
        ([{"=A": 1.0}, {"+G|+G|+G|*AT": 0.9, "-A": 0.1}, {"=T": 1.0}], "AGGGTT"),  # insertion substitution
        ([{"=A": 1.0}, {"+G|+G|+G|N": 0.9, "-A": 0.1}, {"=T": 1.0}], "AGGGNT"),  # insertion N
        ([{"=A": 1.0}], "A"),
    ],
)
def test_call_sequence(cons_percentage_by_key, expected_sequence):
    result_sequence = _call_sequence(cons_percentage_by_key)
    assert result_sequence == expected_sequence
