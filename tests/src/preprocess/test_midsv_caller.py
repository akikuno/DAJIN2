import pytest
import midsv

from DAJIN2.core.preprocess.midsv_caller import _has_inversion_in_splice
from DAJIN2.core.preprocess.midsv_caller import extract_qname_of_map_ont
from DAJIN2.core.preprocess.midsv_caller import replace_internal_n_to_d
from DAJIN2.core.preprocess.midsv_caller import convert_flag_to_strand

from pathlib import Path
from DAJIN2.core.preprocess.mapping import to_sam


###########################################################
# _has_inversion_in_splice
###########################################################
def test_has_inversion_in_splice():
    # Test cases where there is an inversion in splice (insertion followed by deletion)
    assert _has_inversion_in_splice("4M1I4N")
    assert _has_inversion_in_splice("10M1I5N")
    assert _has_inversion_in_splice("1I1N")

    # Test cases where there is no inversion in splice
    assert not _has_inversion_in_splice("10M")
    assert not _has_inversion_in_splice("5N5M")
    assert not _has_inversion_in_splice("1I1M")


def test_has_inversion_in_splice_random_inversion():
    path_reference = Path("tests", "data", "preprocess", "midsv_caller", "reference.fa")
    path_query = Path("tests", "data", "preprocess", "midsv_caller", "query_inversion.fq")
    preset = "splice"
    test = list(to_sam(str(path_reference), str(path_query), preset=preset))
    test = [s.split("\t") for s in test]
    test_cigar = [s[5] for s in test if not s[0].startswith("@")]
    assert _has_inversion_in_splice(test_cigar[0])


def test_has_inversion_in_splice_random_deletion():
    path_reference = Path("tests", "data", "preprocess", "midsv_caller", "reference.fa")
    path_query = Path("tests", "data", "preprocess", "midsv_caller", "query_deletion.fq")
    preset = "splice"
    test = list(to_sam(str(path_reference), str(path_query), preset=preset))
    test = [s.split("\t") for s in test]
    test_cigar = [s[5] for s in test if not s[0].startswith("@")]
    assert not _has_inversion_in_splice(test_cigar[0])


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
#     assert not _has_inversion_in_splice(test_cigar[0])

###########################################################
# extract_qname_of_map_ont
###########################################################


def test_extract_qname_of_map_ont_simulation():
    sam_ont = [["@header"], ["read1", "", "", "", "", "5M"]]
    sam_splice = [["@header"], ["read1", "", "", "", "", "5M"]]
    qname_of_map_ont = extract_qname_of_map_ont(iter(sam_ont), iter(sam_splice))
    assert qname_of_map_ont == {"read1"}

    # Large deletion
    sam_ont = [["@header"], ["read1", "", "", "", "", "5M"], ["read1", "", "", "", "", "5M"]]
    sam_splice = [["@header"], ["read1", "", "", "", "", "5M100D5M"]]
    qname_of_map_ont = extract_qname_of_map_ont(iter(sam_ont), iter(sam_splice))
    assert qname_of_map_ont == set()

    # Inversion
    sam_ont = [
        ["@header"],
        ["read1", "", "", "", "", "5M"],
        ["read1", "", "", "", "", "100M"],
        ["read1", "", "", "", "", "5M"],
    ]
    sam_splice = [["@header"], ["read1", "", "", "", "", "5M100I100N5M"]]
    qname_of_map_ont = extract_qname_of_map_ont(iter(sam_ont), iter(sam_splice))
    assert qname_of_map_ont == {"read1"}

    # Inversion with single read in map-ont
    sam_ont = [["@header"], ["read1", "", "", "", "", "5M10N5M"]]
    sam_splice = [["@header"], ["read1", "", "", "", "", "5M100I100N5M"]]
    qname_of_map_ont = extract_qname_of_map_ont(iter(sam_ont), iter(sam_splice))
    assert qname_of_map_ont == {"read1"}


def test_extract_qname_of_map_ont_real():
    sam_ont = list(midsv.read_sam(Path("tests", "data", "preprocess", "midsv_caller", "stx2-ont_deletion.sam")))
    sam_splice = list(midsv.read_sam(Path("tests", "data", "preprocess", "midsv_caller", "stx2-splice_deletion.sam")))
    qname_of_map_ont = extract_qname_of_map_ont(iter(sam_ont), iter(sam_splice))
    assert qname_of_map_ont == {"stx2-small-deletion"}


###########################################################
# replace n to d
###########################################################


def test_replace_internal_n_to_d():
    sequence = "XYZABCDEF"

    # 既存のテストケース
    midsv_samples = [{"CSSPLIT": "N,N,A,B,N,N"}]
    result = list(replace_internal_n_to_d(midsv_samples, sequence))
    assert result == [{"CSSPLIT": "N,N,A,B,N,N"}]

    midsv_samples = [{"CSSPLIT": "A,B,N,N,C,D"}]
    result = list(replace_internal_n_to_d(midsv_samples, sequence))
    assert result == [{"CSSPLIT": "A,B,-Z,-A,C,D"}]

    # 新しいテストケース
    # 1. Nが連続していない場合の置換
    midsv_samples = [{"CSSPLIT": "A,N,B,N,C"}]
    result = list(replace_internal_n_to_d(midsv_samples, sequence))
    assert result == [{"CSSPLIT": "A,-Y,B,-A,C"}]

    # 2. Nが最初と最後に連続している場合
    midsv_samples = [{"CSSPLIT": "N,N,N,A,B,C,N,N"}]
    result = list(replace_internal_n_to_d(midsv_samples, sequence))
    assert result == [{"CSSPLIT": "N,N,N,A,B,C,N,N"}]

    # 3. Nが複数回連続している場合の置換
    midsv_samples = [{"CSSPLIT": "A,B,N,N,C,N,N,D"}]
    result = list(replace_internal_n_to_d(midsv_samples, sequence))
    assert result == [{"CSSPLIT": "A,B,-Z,-A,C,-C,-D,D"}]

    # 4. Nが1つだけ存在する場合の置換
    midsv_samples = [{"CSSPLIT": "A,B,N,C,D"}]
    result = list(replace_internal_n_to_d(midsv_samples, sequence))
    assert result == [{"CSSPLIT": "A,B,-Z,C,D"}]

    # 5. Nが存在しない場合の確認
    midsv_samples = [{"CSSPLIT": "A,B,C,D,E"}]
    result = list(replace_internal_n_to_d(midsv_samples, sequence))
    assert result == [{"CSSPLIT": "A,B,C,D,E"}]


def test_replace_internal_n_to_d_multiple_samples():
    midsv_sample = [{"CSSPLIT": "N,N,N,=A,N,=C,N,N"}, {"CSSPLIT": "N,N,=A,N,N,=C,=C,N"}]
    sequence = "GCAACCCC"
    test = replace_internal_n_to_d(midsv_sample, sequence)
    test = list(test)
    answer = [{"CSSPLIT": "N,N,N,=A,-C,=C,N,N"}, {"CSSPLIT": "N,N,=A,-A,-C,=C,=C,N"}]
    assert test == answer


def test_replace_internal_n_to_d_large_n():
    midsv_sample = [{"CSSPLIT": "N,N,N,N,N,N,=C,N,=A"}]
    sequence = "GCAACCCCA"
    test = replace_internal_n_to_d(midsv_sample, sequence)
    test = list(test)
    answer = [{"CSSPLIT": "N,N,N,N,N,N,=C,-C,=A"}]
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
