from collections import defaultdict
from typing import NamedTuple

import pytest
from DAJIN2.core.consensus.name_handler import (
    _detect_sv,
    _determine_suffix,
    _format_allele_label,
    add_key_by_allele_name,
    call_allele_name,
    update_key_by_allele_name,
)

###########################################################
# detect_sv
###########################################################


def test_detect_sv_threshold():
    cons_percentages = defaultdict(list)
    cons_percentages["sample1"] = [{"N": 90, "=C": 10}]  # one n
    cons_percentages["sample2"] = [{"+G|+G|=A": 80, "=C": 20}]  # two insertion
    cons_percentages["sample3"] = [{"-A": 100}, {"-A": 100}, {"-A": 100}]  # three deletion
    cons_percentages["sample4"] = [{"*AT": 100}, {"*AT": 100}, {"*AT": 100}, {"*AT": 100}]  # four substitution
    cons_percentages["sample5"] = [{"=a": 100}]  # inversion
    assert _detect_sv(cons_percentages, threshold=1) == [True, True, True, True, True]
    assert _detect_sv(cons_percentages, threshold=2) == [False, True, True, True, True]
    assert _detect_sv(cons_percentages, threshold=3) == [False, False, True, True, True]
    assert _detect_sv(cons_percentages, threshold=4) == [False, False, False, True, True]
    assert _detect_sv(cons_percentages, threshold=5) == [False, False, False, False, True]
    assert _detect_sv(cons_percentages, threshold=6) == [False, False, False, False, True]


###########################################################
# call_allele_name
###########################################################


@pytest.mark.parametrize(
    "label, total_labels, expected_output",
    [
        (1, 10, "01"),
        (5, 100, "005"),
        (10, 10, "10"),
        (99, 99, "99"),
        (1, 1000, "0001"),
    ],
)
def test_format_allele_label(label, total_labels, expected_output):
    result = _format_allele_label(label, total_labels)
    assert result == expected_output


# Test for determine_suffix function
@pytest.mark.parametrize(
    "cons_seq, fasta_allele, is_sv, expected_output",
    [
        ("ATCG", "ATCG", False, "_intact"),
        ("ATCG", "ATCC", True, "_sv"),
        ("ATCG", "ATCC", False, "_indels"),
    ],
)
def test_determine_suffix(cons_seq, fasta_allele, is_sv, expected_output):
    result = _determine_suffix(cons_seq, fasta_allele, is_sv)
    assert result == expected_output


class ConsensusKey(NamedTuple):
    allele: str
    label: int
    percent: float


# Example test cases for call_allele_name function
@pytest.mark.parametrize(
    "cons_sequences, cons_percentages, FASTA_ALLELES, threshold, expected_output",
    [
        # Here, you can add test cases with the corresponding expected output.
        # (cons_sequences, cons_percentages, FASTA_ALLELES, threshold, expected_output)
        (
            {ConsensusKey("control", 1, 100): "ACGT"},
            {ConsensusKey("control", 1, 100): [{"A": 100}, {"C": 100}, {"G": 100}, {"T": 100}]},
            {"control": "ACGT"},
            50,
            {1: "allele1_control_intact_100%"},
        ),
        (
            {ConsensusKey("control", 10, 100): "ACGT"},
            {ConsensusKey("control", 10, 100): [{"A": 100}, {"C": 100}, {"G": 100}, {"T": 100}]},
            {"control": "ACGT"},
            50,
            {10: "allele10_control_intact_100%"},
        ),
    ],
)
def test_call_allele_name(cons_sequences, cons_percentages, FASTA_ALLELES, threshold, expected_output):
    result = call_allele_name(cons_sequences, cons_percentages, FASTA_ALLELES, threshold)
    assert result == expected_output


# Example test cases for cakey_by_allele_name function
@pytest.mark.parametrize(
    "cons, allele_names, expected_output",
    [
        (
            {ConsensusKey("control", 1, 100): "value1", ConsensusKey("control", 2, 100): "value2"},
            {1: "name1", 2: "name2"},
            {"name1": "value1", "name2": "value2"},
        ),
    ],
)
def test_update_key_by_allele_name(cons, allele_names, expected_output):
    result = update_key_by_allele_name(cons, allele_names)
    assert result == expected_output


# Example test cases for add_key_by_allele_name function
@pytest.mark.parametrize(
    "clust_sample, allele_names, expected_output",
    [
        ([{"LABEL": 1}], {1: "name1"}, [{"LABEL": 1, "NAME": "name1"}]),
        # Add more test cases
    ],
)
def test_add_key_by_allele_name(clust_sample, allele_names, expected_output):
    result = add_key_by_allele_name(clust_sample, allele_names)
    assert result == expected_output
