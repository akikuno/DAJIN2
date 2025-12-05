from __future__ import annotations

from pathlib import Path

import pytest

from DAJIN2.core.preprocess.alignment.mapping import to_sam
from DAJIN2.core.preprocess.external_integration.midsv_caller import (
    convert_consecutive_indels_to_match,
    convert_flag_to_strand,
    has_inversion_in_splice,
    replace_internal_n_to_d,
)


###########################################################
# has_inversion_in_splice
###########################################################
def test_has_inversion_in_splice():
    # Test cases where there is an inversion in splice (insertion followed by deletion)
    assert has_inversion_in_splice("4M1I4N")
    assert has_inversion_in_splice("10M1I5N")
    assert has_inversion_in_splice("1I1N")

    # Test cases where there is no inversion in splice
    assert not has_inversion_in_splice("10M")
    assert not has_inversion_in_splice("5N5M")
    assert not has_inversion_in_splice("1I1M")


def test_has_inversion_in_splice_random_inversion():
    path_reference = Path("tests", "data", "preprocess", "midsv_caller", "reference.fa")
    path_query = Path("tests", "data", "preprocess", "midsv_caller", "query_inversion.fq")
    preset = "splice"
    test = list(to_sam(str(path_reference), str(path_query), preset=preset))
    test = [s.split("\t") for s in test]
    test_cigar = [s[5] for s in test if not s[0].startswith("@")]
    assert has_inversion_in_splice(test_cigar[0])


def test_has_inversion_in_splice_random_deletion():
    path_reference = Path("tests", "data", "preprocess", "midsv_caller", "reference.fa")
    path_query = Path("tests", "data", "preprocess", "midsv_caller", "query_deletion.fq")
    preset = "splice"
    test = list(to_sam(str(path_reference), str(path_query), preset=preset))
    test = [s.split("\t") for s in test]
    test_cigar = [s[5] for s in test if not s[0].startswith("@")]
    assert not has_inversion_in_splice(test_cigar[0])


"""
- insertionはmap-ontで１本のリードとして表現されるので、has_inversion_in_spliceは実行されない
"""
# def test_has_inversion_in_splice_random_insertion():
#     path_reference = Path("tests", "data", "preprocess", "midsv_caller", "reference.fa")
#     path_query = Path("tests", "data", "preprocess", "midsv_caller", "query_insertion.fq")
#     preset = "splice"
#     preset = "map-ont"
#     test = list(to_sam(str(path_reference), str(path_query), preset=preset))
#     test = [s.split("\t") for s in test]
#     test_cigar = [s[5] for s in test if not s[0].startswith("@")]
#     assert not has_inversion_in_splice(test_cigar[0])


###########################################################
# replace n to d
###########################################################


def test_replace_internal_n_to_d():
    sequence = "XYZABCDEF"

    # 既存のテストケース
    midsv_samples = [{"MIDSV": "N,N,A,B,N,N"}]
    result = list(replace_internal_n_to_d(midsv_samples, sequence))
    assert result == [{"MIDSV": "N,N,A,B,N,N"}]

    midsv_samples = [{"MIDSV": "A,B,N,N,C,D"}]
    result = list(replace_internal_n_to_d(midsv_samples, sequence))
    assert result == [{"MIDSV": "A,B,-Z,-A,C,D"}]

    # 新しいテストケース
    # 1. Nが連続していない場合の置換
    midsv_samples = [{"MIDSV": "A,N,B,N,C"}]
    result = list(replace_internal_n_to_d(midsv_samples, sequence))
    assert result == [{"MIDSV": "A,-Y,B,-A,C"}]

    # 2. Nが最初と最後に連続している場合
    midsv_samples = [{"MIDSV": "N,N,N,A,B,C,N,N"}]
    result = list(replace_internal_n_to_d(midsv_samples, sequence))
    assert result == [{"MIDSV": "N,N,N,A,B,C,N,N"}]

    # 3. Nが複数回連続している場合の置換
    midsv_samples = [{"MIDSV": "A,B,N,N,C,N,N,D"}]
    result = list(replace_internal_n_to_d(midsv_samples, sequence))
    assert result == [{"MIDSV": "A,B,-Z,-A,C,-C,-D,D"}]

    # 4. Nが1つだけ存在する場合の置換
    midsv_samples = [{"MIDSV": "A,B,N,C,D"}]
    result = list(replace_internal_n_to_d(midsv_samples, sequence))
    assert result == [{"MIDSV": "A,B,-Z,C,D"}]

    # 5. Nが存在しない場合の確認
    midsv_samples = [{"MIDSV": "A,B,C,D,E"}]
    result = list(replace_internal_n_to_d(midsv_samples, sequence))
    assert result == [{"MIDSV": "A,B,C,D,E"}]


def test_replace_internal_n_to_d_multiple_samples():
    midsv_sample = [{"MIDSV": "N,N,N,=A,N,=C,N,N"}, {"MIDSV": "N,N,=A,N,N,=C,=C,N"}]
    sequence = "GCAACCCC"
    test = replace_internal_n_to_d(midsv_sample, sequence)
    test = list(test)
    answer = [{"MIDSV": "N,N,N,=A,-C,=C,N,N"}, {"MIDSV": "N,N,=A,-A,-C,=C,=C,N"}]
    assert test == answer


def test_replace_internal_n_to_d_large_n():
    midsv_sample = [{"MIDSV": "N,N,N,N,N,N,=C,N,=A"}]
    sequence = "GCAACCCCA"
    test = replace_internal_n_to_d(midsv_sample, sequence)
    test = list(test)
    answer = [{"MIDSV": "N,N,N,N,N,N,=C,-C,=A"}]
    assert test == answer


###########################################################
# convert_flag_to_strand
###########################################################
@pytest.mark.parametrize(
    "input_sample, expected_output",
    [
        # FLAGが16の場合
        ([{"FLAG": 16}], [{"STRAND": "-"}]),
        # FLAGが2064の場合
        ([{"FLAG": 2064}], [{"STRAND": "-"}]),
        # FLAGが上記のいずれでもない場合
        ([{"FLAG": 32}], [{"STRAND": "+"}]),
        # 複数のサンプルをテスト
        ([{"FLAG": 16}, {"FLAG": 2064}, {"FLAG": 32}], [{"STRAND": "-"}, {"STRAND": "-"}, {"STRAND": "+"}]),
    ],
)
def test_convert_flag_to_strand(input_sample, expected_output):
    result = list(convert_flag_to_strand(iter(input_sample)))
    assert result == expected_output


###########################################################
# convert_consecutive_indels_to_match
###########################################################


@pytest.mark.parametrize(
    "cons, expected",
    [
        # simple case
        ("=A,-C,-G,-T,+T|=A", "=A,-C,-G,=T,=A"),
        ("=A,-C,-G,-T,+G|+T|=A", "=A,-C,=G,=T,=A"),
        ("=A,-C,-G,-T,+C|+G|+T|=A", "=A,=C,=G,=T,=A"),
        # no change
        ("=A,=C,=G,=T,=A", "=A,=C,=G,=T,=A"),
        # empty
        ("", ""),
    ],
)
def test_convert_consecutive_indels_to_match(cons, expected):
    assert convert_consecutive_indels_to_match(cons) == expected
