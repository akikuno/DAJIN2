import pytest

from src.DAJIN2.core.consensus import sv_annotator


@pytest.mark.parametrize(
    "cons_tag_midsv, midsv_sv_allele, expected",
    [
        # Deletion
        (["=A", "=A"], ["=A", "-T", "-T", "=A"], ["=A", "-T", "-T", "=A"]),
        (
            ["=A", "=A", "=A", "=A"],
            ["=A", "-T", "-T", "=A", "-C", "-C", "=A"],
            ["=A", "-T", "-T", "=A", "-C", "-C", "=A"],
        ),
        # Inversion
        (["=A", "=C", "=C", "=A"], ["=A", "=c", "=c", "=A"], ["=A", "=c", "=c", "=A"]),
        (["=A", "-C", "*CG", "+T|=C", "=A"], ["=A", "=c", "=c", "=c", "=A"], ["=A", "-c", "*cg", "+t|=c", "=A"]),
        # Insertion
        (["=A", "=C", "=C", "=G", "=T"], ["=A", "+C|+C|=G", "=T"], ["=A", "+C|+C|=G", "=T"]),
        (["=A", "=C", "-C", "=G", "=T"], ["=A", "+C|+C|=G", "=T"], ["=A", "+C|=G", "=T"]),
        (["=A", "=C", "=C", "-G", "=T"], ["=A", "+C|+C|=G", "=T"], ["=A", "+C|+C|=T"]),
        (["=A", "=C", "=C", "-G", "-T", "=A"], ["=A", "+C|+C|=G", "=T", "=A"], ["=A", "+C|+C|=A"]),
        # Complex insertion
        (
            ["=A", "=C", "=C", "=G", "=A", "=C", "=C", "=G", "=A"],
            ["=A", "+C|+C|=G", "=A", "+C|+C|=G", "=A"],
            ["=A", "+C|+C|=G", "=A", "+C|+C|=G", "=A"],
        ),
    ],
)
def test_annotate_sv_allele(cons_tag_midsv, midsv_sv_allele, expected):
    assert sv_annotator.annotate_sv_allele(cons_tag_midsv, midsv_sv_allele) == expected


###############################################################################
# reflect_sv_deletion_in_midsv
###############################################################################


@pytest.mark.parametrize(
    "sv_midsv_tag, expected",
    [
        (["=A", "-N", "=T", "=T", "-N", "=C", "=C", "-N", "=A", "=A"], [(2, 3), (5, 6), ()]),
        (["-N", "=G", "=G", "-N", "=T", "-N", "=C"], [(1, 2), (4, 4), ()]),
        (["=A", "=T", "=G"], [()]),  # No deletions, no flanked regions
        ([], [()]),  # Empty input
    ],
)
def test_extract_flanked_ranges(sv_midsv_tag, expected):
    assert sv_annotator.extract_flanked_ranges(sv_midsv_tag) == expected


@pytest.mark.parametrize(
    "tag_flanked, is_insertions, expected",
    [
        # Case 1: Basic conversion
        ([["=A", "-N", "=T"], ["=C", "=G"], ["-N", "=A"]], [True, False, True], [["+A", "+T"], ["=C", "=G"], ["+A"]]),
        # Case 2: No conversion (all False)
        (
            [["=A", "=T"], ["-N", "=C"], ["=G", "-N"]],
            [False, False, False],
            [["=A", "=T"], ["-N", "=C"], ["=G", "-N"]],
        ),
        # Case 3: All elements are insertions
        ([["=A", "=T"], ["=C", "=G"], ["=G", "=A"]], [True, True, True], [["+A", "+T"], ["+C", "+G"], ["+G", "+A"]]),
        # Case 4: List with only deletions
        (
            [["-N", "-N"], ["-X", "-Y"], ["-A", "-B"]],
            [True, True, True],
            [[], [], []],
        ),  # All elements are removed since they are deletions
        # Case 5: Empty list
        ([], [], []),
        # Case 6: Mixed tags with insertions
        (
            [["=A", "-N", "=T"], ["-N", "=G"], ["=C", "=G", "-N"]],
            [False, True, True],
            [["=A", "-N", "=T"], ["+G"], ["+C", "+G"]],
        ),
    ],
)
def test_convert_tags_to_insertion(tag_flanked, is_insertions, expected):
    assert sv_annotator.convert_tags_to_insertion(tag_flanked, is_insertions) == expected


@pytest.mark.parametrize(
    "tag_flanked, is_insertions, expected",
    [
        # Test case 1: Merging consecutive insertions
        ([["+G", "+G"], ["+C"]], [True, True, False], [[], [], ["+G", "+G", "+C"]]),
        # Test case 2: Single insertion region
        ([["+A"]], [True, False, False], [[], ["+A"]]),
        # Test case 3: No insertions, should return None for all
        ([["=T"], ["=C"]], [False, False, False], [None, None]),
        # Test case 4: Mixed insertions and normal regions
        ([["+A"], ["=T"], ["+C"]], [True, False, True, False], [[], ["+A"], None, [], ["+C"]]),
        # Test case 5: Empty input
        ([], [], []),
        # Test case 6: Multiple insertions merging correctly
        ([["+T"], ["+A"], ["+C"], ["=G"]], [True, True, True, False, False], [[], [], [], ["+T", "+A", "+C"], None]),
        # Test case 7: Insertions separated by a normal region
        ([["+T"], ["=A"], ["+C"]], [True, False, True, False], [[], ["+T"], None, [], ["+C"]]),
    ],
)
def test_merge_inserted_sequences(tag_flanked, is_insertions, expected):
    assert sv_annotator.merge_inserted_sequences(tag_flanked, is_insertions) == expected


@pytest.mark.parametrize(
    "sv_midsv_tag, expected",
    [
        # Test case 1: Basic case with multiple deletions
        (["=A", "-N", "-N", "-N", "=T", "=T", "-N", "-N", "=C", "=C", "-N", "=A", "=A"], [4, 8, 11]),
        # Test case 2: Single deletion
        (["-N", "=A"], [1]),
        # Test case 3: No deletions
        (["=A", "=T", "=C"], []),
        # Test case 4: Deletion at the end
        (["=A", "-N", "-N", "=T", "-N"], [3]),
        # Test case 5: Consecutive deletions at different positions
        (["-N", "-N", "=A", "-N", "=C", "-N", "-N", "=G"], [2, 4, 7]),
        # Test case 6: Empty input
        ([], []),
        # Test case 7: Deletions at the start
        (["-N", "-N", "-N", "=A", "=T"], [3]),
    ],
)
def test_get_end_of_deletion_indices(sv_midsv_tag, expected):
    assert sv_annotator.get_end_of_deletion_indices(sv_midsv_tag) == expected


@pytest.mark.parametrize(
    "cons_midsv_tag, sv_midsv_tag, tag_insertions, index_flanked, index_end_of_deletion, expected",
    [
        # Test case 1: Basic case with insertions and deletions
        (
            ["=A", "-N", "*TG", "*TG", "-N", "*TC", "-C", "-N", "=A", "=A"],
            ["=A", "-N", "=T", "=T", "-N", "=T", "=C", "-N", "=A", "=A"],
            [[], [], ["+G", "+G", "+C"]],
            [(2, 3), (5, 6), ()],
            [2, 5, 8],
            ["=A", "-N", "-T", "-T", "-N", "-T", "-C", "-N", "+G|+G|+C|=A", "=A"],
        ),
        # Test case 2: No insertions
        (
            ["=A", "-N", "=T", "=C", "-N", "=G"],
            ["=A", "-N", "=T", "=C", "-N", "=G"],
            [None, None],
            [(2, 3), ()],
            [2, 5],
            ["=A", "-N", "=T", "=C", "-N", "=G"],
        ),
        # Test case 3: Single insertion region
        (
            ["=A", "-N", "=G", "-N", "=T"],
            ["=A", "-N", "=G", "-N", "=T"],
            [["+C"], None, None],
            [[0, 1], (), ()],
            [0, 3, 4],
            ["+C|=A", "-N", "=G", "-N", "=T"],
        ),
        # Test case 4: Empty list
        ((), (), (), (), (), ()),
    ],
)
def test_annotate_insertion(
    cons_midsv_tag, sv_midsv_tag, tag_insertions, index_flanked, index_end_of_deletion, expected
):
    assert (
        sv_annotator.annotate_insertion(
            cons_midsv_tag, sv_midsv_tag, tag_insertions, index_flanked, index_end_of_deletion
        )
        == expected
    )
