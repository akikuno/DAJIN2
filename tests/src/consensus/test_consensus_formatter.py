import pytest

from DAJIN2.core.consensus.consensus_formatter import (
    call_allele_name,
    merge_duplicated_cons_others,
    merge_duplicated_cons_sequences,
    update_key_by_allele_name,
    update_label_percent_readnum_name,
)
from DAJIN2.utils.config import ConsensusKey

###########################################################
# call_allele_name
###########################################################


@pytest.mark.parametrize(
    "cons_sequences, fasta_alleles, expected_label_name, expected_name_allele",
    [
        (
            {
                ConsensusKey("control", 1, 75): "ATCG",
                ConsensusKey("mutant", 2, 25): "ATCC",
            },
            {"control": "ATCG", "mutant": "ATCG"},
            {
                1: "allele01|control|intact|75%",
                2: "allele02|mutant|indels|25%",
            },
            {
                "allele01|control|intact|75%": "control",
                "allele02|mutant|indels|25%": "mutant",
            },
        ),
        (
            {
                ConsensusKey("control", 1, 90): "ATCG",
                ConsensusKey("deletion3_DAJIN2predicted", 3, 10): "ATCG",
            },
            {"control": "ATCG", "deletion3_DAJIN2predicted": "ATCG"},
            {
                1: "allele01|control|intact|90%",
                3: "allele03|unassigned|deletion|10%",
            },
            {
                "allele01|control|intact|90%": "control",
                "allele03|unassigned|deletion|10%": "deletion3_DAJIN2predicted",
            },
        ),
    ],
)
def test_call_allele_name(cons_sequences, fasta_alleles, expected_label_name, expected_name_allele):
    got_label_name, got_name_allele = call_allele_name(cons_sequences, fasta_alleles)
    assert got_label_name == expected_label_name
    assert got_name_allele == expected_name_allele


def test_call_allele_name_digit_width():
    cons_sequences = {ConsensusKey("control", i, 1): "ATCG" for i in range(1, 101)}
    fasta_alleles = {"control": "ATCG"}

    got_label_name, got_name_allele = call_allele_name(cons_sequences, fasta_alleles)

    assert got_label_name[1] == "allele001|control|intact|1%"
    assert got_label_name[100] == "allele100|control|intact|1%"
    assert got_name_allele["allele001|control|intact|1%"] == "control"


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
            {1: "name1|30%", 3: "name3|10%"},
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
