from typing import NamedTuple

import pytest

from DAJIN2.core.consensus.consensus_formatter import (
    ConsensusKey as ConsensusKeyData,
    call_allele_name,
    detect_sv,
    determine_suffix,
    format_allele_label,
    generate_allele_mapping,
    merge_duplicated_cons_others,
    merge_duplicated_cons_sequences,
    update_key_by_allele_name,
    update_label_percent_readnum_name,
)

###########################################################
# detect_sv
###########################################################


@pytest.mark.parametrize(
    "cons_per, threshold, expected",
    [
        ([{"=N": 90, "=C": 10}], 1, True),  # one "=N"
        ([{"+G|+G|=A": 80, "=C": 20}], 1, True),  # two insertions
        ([{"-A": 100}, {"-A": 100}, {"-A": 100}], 3, True),  # three deletions
        ([{"*AT": 100}, {"*AT": 100}, {"*AT": 100}, {"*AT": 100}], 4, True),  # four substitutions
        ([{"=a": 100}], 5, True),  # inversion
        ([{"=N": 90, "=C": 10}], 2, False),  # fails threshold
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


def test_call_allele_name_uses_sv_display_map():
    key_sv = ConsensusKeyData(allele="deletion001__uuid", label=1, percent=60.0)
    cons_sequences = {key_sv: "A"}
    cons_percentages = {key_sv: [{"-A": 100.0}]}
    fasta_alleles = {"deletion001__uuid": "A"}
    sv_name_map = {"deletion001__uuid": "DAJIN_deletion01"}

    allele_names, final_sv_map = call_allele_name(
        cons_sequences, cons_percentages, fasta_alleles, sv_threshold=1, sv_name_map=sv_name_map
    )

    assert allele_names[1].startswith("allele01_unintended_deletion_1_60")
    assert final_sv_map == {"deletion001__uuid": "unintended_deletion_1"}


def test_call_allele_name_non_sv_readable_labels():
    key_control = ConsensusKeyData(allele="control", label=1, percent=70.0)
    key_sample = ConsensusKeyData(allele="flox", label=2, percent=30.0)
    cons_sequences = {key_control: "AC", key_sample: "AT"}
    cons_percentages = {key_control: [{"=A": 100.0}], key_sample: [{"*AC": 100.0}]}
    fasta_alleles = {"control": "AC", "flox": "AC"}

    allele_names, _ = call_allele_name(cons_sequences, cons_percentages, fasta_alleles, sv_threshold=1)

    assert allele_names[1].startswith("allele01_control_70")
    assert allele_names[2].startswith("allele02_flox_with_indels_30")


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
        # Case 6: Consecutive alleles
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
    "clust_sample, allele_names, label_before_to_after, expected_output",
    [
        (
            [
                {"LABEL": 1, "PERCENT": 20, "READNUM": 30},
                {"LABEL": 2, "PERCENT": 10, "READNUM": 20},
                {"LABEL": 3, "PERCENT": 10, "READNUM": 10},
            ],
            {1: "name1_30%", 3: "name3_10%"},
            {2: 1},
            [
                {"LABEL": 1, "NAME": "name1_30%", "PERCENT": 30.0, "READNUM": 50},
                {"LABEL": 1, "NAME": "name1_30%", "PERCENT": 30.0, "READNUM": 50},
                {"LABEL": 3, "NAME": "name3_10%", "PERCENT": 10.0, "READNUM": 10},
            ],
        ),
        # Add more test cases
    ],
)
def test_update_label_percent_readnum_name(clust_sample, allele_names, label_before_to_after, expected_output):
    result = update_label_percent_readnum_name(clust_sample, allele_names, label_before_to_after)
    assert result == expected_output


###########################################################
# Merge alleles where cons_sequences and cons_percentages are an exact match
###########################################################


# === Helpers (same style as provided) ===
def to_tuple_set_str(d: dict[ConsensusKey, str]) -> set[tuple[str, int, float, str]]:
    return {(k.allele, k.label, float(k.percent), v) for k, v in d.items()}


def to_tuple_set_any(d: dict[ConsensusKey, object]) -> set[tuple[str, int, float, object]]:
    return {(k.allele, k.label, float(k.percent), v) for k, v in d.items()}


def dict_of_sets_equal(a: dict[int, set[int]], b: dict[int, set[int]]) -> bool:
    if set(a.keys()) != set(b.keys()):
        return False
    return all(set(a[k]) == set(b[k]) for k in a.keys())


###########################################################
# Merge alleles where cons_sequences and cons_percentages are an exact match
# Then re-rank merged representatives globally by total percent (desc) to labels 1..N.
###########################################################
@pytest.mark.parametrize(
    "inp, expected_merged, expected_map",
    [
        # Case 1:
        # Groups:
        #  - control/AAAA: labels {1(15),3(5)} -> merged total 20, interim=1
        #  - control/BBBB: {2(10)}            -> total 10, interim=2
        #  - case/AAAA:    {4(7)}             -> total 7,  interim=4
        # Ranking by total percent: 20 -> rank1, 10 -> rank2, 7 -> rank3
        # Final merged keys use rank labels; mapping is original -> final rank.
        (
            {
                ConsensusKey(allele="control", label=1, percent=15): "AAAA",
                ConsensusKey(allele="control", label=2, percent=10): "BBBB",
                ConsensusKey(allele="control", label=3, percent=5): "AAAA",
                ConsensusKey(allele="case", label=4, percent=7): "AAAA",
            },
            {
                ConsensusKey(allele="control", label=1, percent=20): "AAAA",  # rank1
                ConsensusKey(allele="control", label=2, percent=10): "BBBB",  # rank2
                ConsensusKey(allele="case", label=3, percent=7): "AAAA",  # rank3
            },
            {1: 1, 3: 1, 2: 2, 4: 3},  # original -> final rank
        ),
        # Case 2:
        # Groups:
        #  - control/AAAA: {1(15),3(5)} -> total 20, interim=1
        #  - control/BBBB: {2(10)}      -> total 10, interim=2
        #  - case/AAAA:    {5(7),10(7),20(7)} -> total 21, interim=5
        # Ranking: 21 -> rank1 (case/AAAA), 20 -> rank2 (control/AAAA), 10 -> rank3 (control/BBBB)
        (
            {
                ConsensusKey(allele="control", label=1, percent=15): "AAAA",
                ConsensusKey(allele="control", label=2, percent=10): "BBBB",
                ConsensusKey(allele="control", label=3, percent=5): "AAAA",
                ConsensusKey(allele="case", label=5, percent=7): "AAAA",
                ConsensusKey(allele="case", label=10, percent=7): "AAAA",
                ConsensusKey(allele="case", label=20, percent=7): "AAAA",
            },
            {
                ConsensusKey(allele="case", label=1, percent=21): "AAAA",  # rank1
                ConsensusKey(allele="control", label=2, percent=20): "AAAA",  # rank2
                ConsensusKey(allele="control", label=3, percent=10): "BBBB",  # rank3
            },
            {1: 2, 3: 2, 2: 3, 5: 1, 10: 1, 20: 1},  # original -> final rank
        ),
        # Case 3: empty input
        ({}, {}, {}),
        # Case 4: no merges (all (allele, seq) unique)
        # Totals: control/AAAA=10, control/BBBB=20, case/CCCC=30
        # Ranking: 30->rank1 (case), 20->rank2 (control/BBBB), 10->rank3 (control/AAAA)
        (
            {
                ConsensusKey(allele="control", label=1, percent=10): "AAAA",
                ConsensusKey(allele="control", label=2, percent=20): "BBBB",
                ConsensusKey(allele="case", label=3, percent=30): "CCCC",
            },
            {
                ConsensusKey(allele="case", label=1, percent=30): "CCCC",  # rank1
                ConsensusKey(allele="control", label=2, percent=20): "BBBB",  # rank2
                ConsensusKey(allele="control", label=3, percent=10): "AAAA",  # rank3
            },
            {1: 3, 2: 2, 3: 1},
        ),
    ],
)
def test_merge_duplicated_cons_sequences_with_ranking(inp, expected_merged, expected_map):
    got_merged, got_map = merge_duplicated_cons_sequences(inp)
    assert to_tuple_set_str(got_merged) == to_tuple_set_str(expected_merged)
    assert got_map == expected_map


###########################################################
# NEW: cons_others update that also applies label_before_to_after
# Behavior:
#  1) Relabel each key in cons_others using label_before_to_after (original -> final rank).
#  2) Drop entries whose (final) label does not exist in cons_sequences_merged.
#  3) Update percent in keys to the percent found in cons_sequences_merged for that final label.
###########################################################
@pytest.mark.parametrize(
    "cons_others, cons_sequences_merged, label_before_to_after, expected",
    [
        # Case A (basic): relabel (1->2, 2->3, 3->1), then align with merged labels {1,2}
        # cons_others:
        #   control/1 -> relabeled to control/2 -> kept; percent becomes 10 (from label 2)
        #   control/2 -> relabeled to control/3 -> dropped (label 3 not in cons_sequences_merged)
        #   case/3    -> relabeled to case/1    -> kept; percent becomes 20 (from label 1)
        (
            {
                ConsensusKey("control", 1, 0): "metaA",
                ConsensusKey("control", 2, 0): "metaA",
                ConsensusKey("case", 3, 0): "metaC",
            },
            {
                ConsensusKey("control", 1, 20): "AAAA",  # final label 1
                ConsensusKey("case", 2, 10): "BBBB",  # final label 2
            },
            {1: 1, 2: 1, 3: 2},
            {
                ConsensusKey("control", 1, 20): "metaA",  # from original 1 -> final 2
                ConsensusKey("case", 2, 10): "metaC",  # from original 3 -> final 1
            },
        ),
        # Case B: cons_others empty -> stays empty regardless of mapping
        (
            {},
            {ConsensusKey("control", 1, 5): "AAAA"},
            {1: 1},
            {},
        ),
        # Case C: after relabel, no labels exist in cons_sequences_merged -> all dropped
        (
            {ConsensusKey("control", 9, 0): "metaX"},
            {ConsensusKey("control", 1, 12): "ZZZZ"},
            {9: 99},  # maps to a label (99) that doesn't exist in merged
            {},
        ),
        # Case D: identity mapping; behaves like the old behavior but with explicit map
        (
            {
                ConsensusKey("control", 1, 0): "metaA",
                ConsensusKey("control", 2, 0): "metaB",
            },
            {
                ConsensusKey("control", 1, 33): "AAAA",
                ConsensusKey("control", 2, 67): "BBBB",
            },
            {1: 1, 2: 2},
            {
                ConsensusKey("control", 1, 33): "metaA",
                ConsensusKey("control", 2, 67): "metaB",
            },
        ),
        # Case E: mixed alleles; ensure allele field is preserved while labels are remapped
        # Here, only the labels present in merged (final=1,2) are kept and updated.
        (
            {
                ConsensusKey("control", 10, 0): "metaCtl",
                ConsensusKey("case", 20, 0): "metaCase",
                ConsensusKey("case", 21, 0): "metaDrop",
            },
            {
                ConsensusKey("control", 1, 40): "SEQ1",
                ConsensusKey("case", 2, 60): "SEQ2",
            },
            {10: 1, 20: 2, 21: 99},
            {
                ConsensusKey("control", 1, 40): "metaCtl",  # control keeps allele, label->1
                ConsensusKey("case", 2, 60): "metaCase",  # case keeps allele, label->2
            },
        ),
    ],
)
def test_merge_duplicated_cons_others(cons_others, cons_sequences_merged, label_before_to_after, expected):
    got = merge_duplicated_cons_others(cons_others, cons_sequences_merged, label_before_to_after)
    assert to_tuple_set_any(got) == to_tuple_set_any(expected)
