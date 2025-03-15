from typing import NamedTuple

import pytest
from DAJIN2.core.consensus.name_handler import (
    add_key_by_allele_name,
    detect_sv,
    determine_suffix,
    format_allele_label,
    generate_allele_mapping,
    update_key_by_allele_name,
)

###########################################################
# detect_sv
###########################################################


@pytest.mark.parametrize(
    "cons_per, threshold, expected",
    [
        ([{"N": 90, "=C": 10}], 1, True),  # one "N"
        ([{"+G|+G|=A": 80, "=C": 20}], 1, True),  # two insertions
        ([{"-A": 100}, {"-A": 100}, {"-A": 100}], 3, True),  # three deletions
        ([{"*AT": 100}, {"*AT": 100}, {"*AT": 100}, {"*AT": 100}], 4, True),  # four substitutions
        ([{"=a": 100}], 5, True),  # inversion
        ([{"N": 90, "=C": 10}], 2, False),  # fails threshold
    ],
)
def test_detect_sv(cons_per, threshold, expected):
    assert detect_sv(cons_per, threshold) == expected


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
    result = format_allele_label(label, total_labels)
    assert result == expected_output


# Test for determine_suffix function
@pytest.mark.parametrize(
    "cons_seq, fasta_allele, is_sv, expected_output",
    [
        ("ATCG", "ATCG", False, "intact"),
        ("ATCG", "ATCC", True, "SV"),
        ("ATCG", "ATCC", False, "indels"),
    ],
)
def test_determine_suffix(cons_seq, fasta_allele, is_sv, expected_output):
    result = determine_suffix(cons_seq, fasta_allele, is_sv)
    assert result == expected_output


@pytest.mark.parametrize(
    "alleles, expected",
    [
        # Case 1: Standard mapping
        (
            ["deletion02", "control", "deletion04", "inversion05", "insertion11"],
            {
                "deletion02": "deletion01",
                "deletion04": "deletion02",
                "inversion05": "inversion01",
                "insertion11": "insertion01",
            },
        ),
        # Case 2: Only one allele in each group
        (
            ["deletion01", "inversion01", "insertion01"],
            {"deletion01": "deletion01", "inversion01": "inversion01", "insertion01": "insertion01"},
        ),
        # Case 3: No valid alleles
        (["control", "unknown"], {}),
        # Case 4: Multiple alleles in a single group
        (
            ["deletion01", "deletion02", "deletion03"],
            {"deletion01": "deletion01", "deletion02": "deletion02", "deletion03": "deletion03"},
        ),
        # Case 5: Mixed valid and invalid alleles
        (
            ["deletion02", "control", "inversion05", "insertion11", "unknown"],
            {"deletion02": "deletion01", "inversion05": "inversion01", "insertion11": "insertion01"},
        ),
        # Case 6: Conseqtive alleles
        (["deletion01", "deletion01", "deletion04"], {"deletion01": "deletion01", "deletion04": "deletion02"}),
    ],
)
def test_generate_allele_mapping(alleles, expected):
    assert generate_allele_mapping(alleles) == expected


class ConsensusKey(NamedTuple):
    allele: str
    label: int
    percent: float


###########################################################
# Replace new allele names to the consensus dictionary
###########################################################


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


###########################################################
# Add `NAME` key to RESULT_SAMPLE
###########################################################


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
