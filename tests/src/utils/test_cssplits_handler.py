from __future__ import annotations

import pytest

from DAJIN2.utils import cssplits_handler

###########################################################
# find_n_boundaries
###########################################################


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
    assert cssplits_handler.find_n_boundaries(cssplits) == expected


###########################################################
# convert cssplits to DNA sequence
###########################################################


@pytest.mark.parametrize(
    "sequence, expected",
    [
        ("AGGCGAacgAaccN", "AGGCGAcgtAggtN"),  # Lowercase inversion
        ("AAAAaaaaAAAA", "AAAAttttAAAA"),  # All lowercase in the middle
        ("aaaa", "tttt"),  # Only lowercase
        ("AGGCGATACC", "AGGCGATACC"),  # No lowercase
        ("a", "t"),  # Single lowercase
        ("", ""),  # Empty string
    ],
)
def test_revcomp_inversion(sequence, expected):
    assert cssplits_handler._revcomp_inversion(sequence) == expected


@pytest.mark.parametrize(
    "cssplits, expected",
    [
        (
            ["=A", "=n", "+a|+t|+g|=a", "=c", "=g", "=A", "=a", "=c", "=c", "=N"],
            "AcgtcatnAggtN",
        ),  # Complex case with insertion and inversion
        (["=A", "=T", "=C", "=G"], "ATCG"),  # Simple match only
        (["+a|+t|+g|=a"], "tcat"),  # Single insertion with match
        (["=a", "=c", "=t"], "agt"),  # Simple inversion
        ([], ""),  # Empty cssplits
    ],
)
def test_convert_cssplits_to_sequence(cssplits, expected):
    assert cssplits_handler.convert_cssplits_to_sequence(cssplits) == expected


###########################################################
# convert cssplits to cstag
###########################################################


@pytest.mark.parametrize(
    "cssplits, expected",
    [
        ([], []),
        (["=A", "=C", "=G"], ["=A", "=C", "=G"]),
        (["N", "=C", "N"], ["=N", "=C", "=N"]),
        (["=A", "+G|+G|+G|N", "=T"], ["=A", "+G|+G|+G|=N", "=T"]),
    ],
)
def test_add_match_operator_to_n(cssplits, expected):
    assert cssplits_handler._add_match_operator_to_n(cssplits) == expected


@pytest.mark.parametrize(
    "input_cssplits, expected",
    [
        ([], []),
        (["=ACGT"], ["=ACGT"]),
        (["=ACGT", "*GC", "=ACGT"], ["=ACGT", "*gc", "=ACGT"]),
    ],
)
def test_standardize_case(input_cssplits, expected):
    assert cssplits_handler._standardize_case(input_cssplits) == expected


@pytest.mark.parametrize(
    "cssplits, expected",
    [
        ([], ""),
        (["=A", "=C", "=G", "=T"], "=ACGT"),
        (["=A", "*GT", "=A"], "=A*gt=A"),
        (["-C", "-G"], "-cg"),
        (["*GC", "*TA"], "*gc*ta"),
        (["+A|+A|=C", "=G"], "+aa=CG"),
        (["+A|+C|=T", "+A|+C|=G"], "+ac=T+ac=G"),
        (["=A", "=C", "+A|+C|=G", "-C", "-G", "*CG"], "=AC+ac=G-cg*cg"),
    ],
)
def test_convert_cssplits_to_cstag(cssplits, expected):
    assert cssplits_handler.convert_cssplits_to_cstag(cssplits) == expected


###########################################################
# call sequence
###########################################################


@pytest.mark.parametrize(
    "cons_percentage, expected_sequence",
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
def test_call_sequence(cons_percentage, expected_sequence):
    assert cssplits_handler.call_sequence(cons_percentage) == expected_sequence



###########################################################
# reflect_sv_deletion_in_midsv
###########################################################

# ---------------------------------------------------------
# highlight_sv_deletions
# ---------------------------------------------------------

@pytest.mark.parametrize(
    "midsv_sv_allele, expected_output",
    [
        # Test case 1: Insertion
        (["=A", "+A|+A|=C", "=G"], ['=A', '+A|+A|=C', '=G']),
        # Test case 2: Consecutive Deletions
        (["=A", "-A", "-G", "=G"], ["=A", "!START_OF_DEL_ALLELE!", "-A", "-G", "!END_OF_DEL_ALLELE!", "=G"]),
        # Test case 3: Single Deletion
        (["=A", "-G", "=G"], ["=A", "!START_OF_DEL_ALLELE!", "-G", "!END_OF_DEL_ALLELE!", "=G"]),
    ],
)
def test_highlight_sv_deletions(midsv_sv_allele, expected_output):
    assert cssplits_handler.highlight_sv_deletions(midsv_sv_allele) == expected_output

# ---------------------------------------------------------
# embed_sv_deletions
# ---------------------------------------------------------

test_data = [
    (["=A", "-T", "-G", "=C", "=T"], ["=A", "=C", "=T"], ['=A', '!START_OF_DEL_ALLELE!', '!END_OF_DEL_ALLELE!', '=C', '=T']),
    (["=A", "=T", "=C", "=T"], ["=A", "-G", "-C", "=T"], ["=A", "-G", "-C", "=T"]),
]
@pytest.mark.parametrize("midsv_sv_allele, midsv_consensus, expected_output", test_data)
def test_embed_sv_deletions(midsv_sv_allele, midsv_consensus, expected_output):
    midsv_sv_allele_highlighted = cssplits_handler.highlight_sv_deletions(midsv_sv_allele)
    result = cssplits_handler.embed_sv_deletions(midsv_consensus, midsv_sv_allele_highlighted)
    assert result == expected_output


# ---------------------------------------------------------
# get_flanked_tags_by_deletions
# ---------------------------------------------------------

@pytest.mark.parametrize(
    "midsv_tags, expected_output",
    [
        # No flanked tags between deletions
        (['=A', '!START_OF_DEL_ALLELE!', '=C', '!END_OF_DEL_ALLELE!', '=T'], ([], [])),

        # One flanked tag between two deletions
        (
            ['=A', '!START_OF_DEL_ALLELE!', '=C', '!END_OF_DEL_ALLELE!', '=T', '!START_OF_DEL_ALLELE!', '=C', '!END_OF_DEL_ALLELE!'],
            ([['=T']], [(4, 5)])
        ),

        # Multiple flanked tags between deletions
        (
            ['=A', '!START_OF_DEL_ALLELE!', '=C', '!END_OF_DEL_ALLELE!', '=T', '!START_OF_DEL_ALLELE!', '=C', '!END_OF_DEL_ALLELE!', '=T', '=T', '!START_OF_DEL_ALLELE!', '=C', '!END_OF_DEL_ALLELE!'],
            ([['=T'], ['=T', '=T']], [(4, 5), (8, 10)])
        ),
    ]
)
def test_get_flanked_tags_by_deletions(midsv_tags, expected_output):
    assert cssplits_handler.get_flanked_tags_by_deletions(midsv_tags) == expected_output


# ---------------------------------------------------------
# has_consecutive_matches
# ---------------------------------------------------------

@pytest.mark.parametrize(
    "flanked_tag, n, expected_output",
    [
        # Case 1: No matching tags
        (["-A", "+T", "-G"], 10, False),

        # Case 2: Less than n consecutive matches
        (["=A", "=T", "=G", "-C", "=A", "=T"], 4, False),

        # Case 3: Exactly n consecutive matches
        (["=A", "=T", "=G", "=C", "=A", "=T", "=G", "=C", "=A", "=T"], 10, True),

        # Case 4: More than n consecutive matches
        (["=A"] * 15, 10, True),

        # Case 5: Consecutive matches but broken before reaching n
        (["=A", "=T", "=G", "=C", "=A", "-T", "=G", "=C", "=A", "=T"], 5, True),

        # Case 6: Edge case where n is 1 (any match should return True)
        (["=A", "-T", "+G"], 1, True),
    ]
)
def test_has_consecutive_matches(flanked_tag, n, expected_output):
    assert cssplits_handler.has_consecutive_matches(flanked_tag, n) == expected_output


# ---------------------------------------------------------
# has_consecutive_matches
# ---------------------------------------------------------

@pytest.mark.parametrize(
    "midsv_tags, expected_output",
    [
        pytest.param(
            ["=A", "!START_OF_DEL_ALLELE!", "=C", "!END_OF_DEL_ALLELE!", "=T"],
            4,
            id="Single deletion block"
        ),
        pytest.param(
            ["=A", "!START_OF_DEL_ALLELE!", "=C", "!END_OF_DEL_ALLELE!", "=T", "!START_OF_DEL_ALLELE!", "=G", "!END_OF_DEL_ALLELE!"],
            8,
            id="Multiple deletion blocks"
        ),
        pytest.param(
            ["=A", "=T", "=G"],
            None,  # Expecting an exception
            id="No deletion block"
        ),
        pytest.param(
            ["=A", "!START_OF_DEL_ALLELE!", "=C"],
            None,  # Expecting an exception
            id="Unclosed deletion block"
        ),
        pytest.param(
            ["!END_OF_DEL_ALLELE!", "=A", "=T"],
            1,
            id="END_OF_DEL_ALLELE at index 0"
        ),
    ]
)
def test_get_last_index_of_del_allele(midsv_tags, expected_output):
    if expected_output is None:
        with pytest.raises(StopIteration):
            cssplits_handler.get_last_index_of_del_allele(midsv_tags)
    else:
        assert cssplits_handler.get_last_index_of_del_allele(midsv_tags) == expected_output, (
            f"Test failed for input: {midsv_tags}\n"
            f"Expected: {expected_output}, but got: {cssplits_handler.get_last_index_of_del_allele(midsv_tags)}"
        )



# ---------------------------------------------------------
# convert_midsv_tags_to_insertion
# ---------------------------------------------------------

@pytest.mark.parametrize(
    "midsv_tags, expected_output",
    [
        pytest.param(
            ["=T", "=A", "+C|+C|=A", "-T", "=G", "=C"],
            "+T|+A|+C|+C|+A|+G|+C",
            id="Standard case with insertions and deletions"
        ),
        pytest.param(
            ["=A", "=C", "=G"],
            "+A|+C|+G",
            id="Only matches"
        ),
        pytest.param(
            ["+A|+T|=C", "=G", "-T", "+A|+G|=C"],
            "+A|+T|+C|+G|+A|+G|+C",
            id="Mixed insertions and deletions"
        ),
        pytest.param(
            ["-A", "-C", "-G"],
            "",
            id="Only deletions"
        ),
        pytest.param(
            ["+A|+C|+G", "+T|+A|=C"],
            "+A|+C|+G|+T|+A|+C",
            id="Only insertions"
        ),
        pytest.param(
            ["=A", "=T", "=G", "-C", "-A", "=C"],
            "+A|+T|+G|+C",
            id="Trailing deletion"
        ),
        pytest.param(
            ["+A|+T|-C", "=G", "=C"],
            "+A|+T|+G|+C",
            id="Insertion ending with deletion"
        ),
        pytest.param(
            [],
            "",
            id="Empty list"
        ),
    ]
)
def test_convert_midsv_tags_to_insertion(midsv_tags, expected_output):
    actual_output = cssplits_handler.convert_midsv_tags_to_insertion(midsv_tags)
    assert actual_output == expected_output, (
        f"Test failed for input: {midsv_tags}\n"
        f"Expected: {expected_output}, but got: {actual_output}"
    )


# ---------------------------------------------------------
# extract_deletion_tags
# ---------------------------------------------------------

@pytest.mark.parametrize(
    "midsv_sv_deletions_highlighted, expected_output",
    [
        pytest.param(
            ["=A", "!START_OF_DEL_ALLELE!", "-C", "-G", "!END_OF_DEL_ALLELE!", "=T"],
            [["-C", "-G"]],
            id="Single deletion block"
        ),
        pytest.param(
            ["=A", "!START_OF_DEL_ALLELE!", "-C", "-G", "!END_OF_DEL_ALLELE!", "=T", "!START_OF_DEL_ALLELE!", "-A", "-T", "!END_OF_DEL_ALLELE!"],
            [["-C", "-G"], ["-A", "-T"]],
            id="Multiple deletion blocks"
        ),
        pytest.param(
            ["=A", "=T", "=G"],
            [],
            id="No deletion blocks"
        ),
        pytest.param(
            ["!START_OF_DEL_ALLELE!", "-C", "-G", "!END_OF_DEL_ALLELE!", "!START_OF_DEL_ALLELE!", "-A", "-T", "!END_OF_DEL_ALLELE!"],
            [["-C", "-G"], ["-A", "-T"]],
            id="Only deletions without matches"
        ),
        pytest.param(
            ["!START_OF_DEL_ALLELE!", "!START_OF_DEL_ALLELE!", "-C", "!END_OF_DEL_ALLELE!"],
            [[], ["-C"]],
            id="Nested deletion start tags"
        ),
        pytest.param(
            [],
            [],
            id="Empty input"
        ),
    ]
)
def test_extract_deletion_tags(midsv_sv_deletions_highlighted, expected_output):
    actual_output = cssplits_handler.extract_deletion_tags(midsv_sv_deletions_highlighted)
    assert actual_output == expected_output, (
        f"Test failed for input: {midsv_sv_deletions_highlighted}\n"
        f"Expected: {expected_output}, but got: {actual_output}"
    )

# ---------------------------------------------------------
# handle_insertions_within_deletions
# ---------------------------------------------------------

@pytest.mark.parametrize(
    "midsv_consensus_with_deletion, deletion_tags, expected_output",
    [
        pytest.param(
            ["=A", "!START_OF_DEL_ALLELE!", "!END_OF_DEL_ALLELE!", "=A", "!START_OF_DEL_ALLELE!", "!END_OF_DEL_ALLELE!", "=A"],
            [["-G", "-G"], ["-T", "-T"]],
            ['=A', '-G', '-G', '-T', '-T', '+A|=A'],
            id="Single deletion block without insertions"
        )
    ]
)
def test_handle_insertions_within_deletions(midsv_consensus_with_deletion, deletion_tags, expected_output):
    flanked_tags, flanked_indices = cssplits_handler.get_flanked_tags_by_deletions(midsv_consensus_with_deletion)
    actual_output = cssplits_handler.handle_insertions_within_deletions(midsv_consensus_with_deletion, flanked_tags, flanked_indices, deletion_tags)
    assert actual_output == expected_output

# ---------------------------------------------------------
# extract_deletion_tags
# ---------------------------------------------------------


@pytest.mark.parametrize(
    "midsv_consensus_with_deletion, deletion_tags, expected_output",
    [
        pytest.param(
            ["=A", "!START_OF_DEL_ALLELE!", "!END_OF_DEL_ALLELE!", "=C","=T"],
            [["-T", "-T"]],
            ["=A", "-T", "-T","=C", "=T"],
            id="Single deletion block without insertions"
        ),
        pytest.param(
            ["=A", "!START_OF_DEL_ALLELE!", "!END_OF_DEL_ALLELE!", "=C", "=G", "=T"],
            [["-T", "-T"]],
            ["=A", "-T", "-T", "=C", "=G","=T"],
            id="Single deletion block with multiple deletions"
        ),
        pytest.param(
            ["=A", "!START_OF_DEL_ALLELE!", "!END_OF_DEL_ALLELE!", "=A", "=A", "!START_OF_DEL_ALLELE!", "!END_OF_DEL_ALLELE!", "=C", "=C"],
            [["-G", "-G", "-G"], ["-T", "-T"]],
            ["=A", "-G", "-G", "-G", "=A", "=A", "-T", "-T", "=C", "=C"],
            id="Multiple deletion blocks"
        ),
    ]
)
def test_reflect_deletions(midsv_consensus_with_deletion, deletion_tags, expected_output):
    actual_output = cssplits_handler.reflect_deletions(midsv_consensus_with_deletion, deletion_tags)
    assert actual_output == expected_output
# ---------------------------------------------------------
# reflect_sv_deletion_in_midsv
# ---------------------------------------------------------
# @pytest.mark.parametrize(
#     "cssplits, expected",
#     [
#         (
#             ["=T"] * 100 + ["-A"] * 300 + ["*TA"] * 10 + ["-A"] * 300 + ["=T"] * 100,
#             ["=T"] * 100
#             + ["-A"] * 300
#             + ["-T"] * 10
#             + ["-A"] * 300
#             + ["+A|+A|+A|+A|+A|+A|+A|+A|+A|+A|=T"]
#             + ["=T"] * 99,
#         ),
#         (
#             ["=T"] * 100 + ["-A"] * 150 + ["=T"] * 10 + ["-A"] * 150 + ["=T"] * 100,
#             ["=T"] * 100 + ["-A"] * 150 + ["=T"] * 10 + ["-A"] * 150 + ["=T"] * 100,
#         ),
#         (
#             ["=T"] * 100
#             + ["-A"] * 100
#             + ["*TA"] * 10
#             + ["-A"] * 100
#             + ["=T"] * 10
#             + ["-A"] * 100
#             + ["*TA"] * 10
#             + ["-A"] * 100
#             + ["=T"] * 100,
#             ["=T"] * 100
#             + ["-A"] * 100
#             + ["-T"] * 10
#             + ["-A"] * 100
#             + ["+A|+A|+A|+A|+A|+A|+A|+A|+A|+A|=T"]
#             + ["=T"] * 9
#             + ["-A"] * 100
#             + ["-T"] * 10
#             + ["-A"] * 100
#             + ["+A|+A|+A|+A|+A|+A|+A|+A|+A|+A|=T"]
#             + ["=T"] * 99,
#         ),
#     ],
#     ids=[
#         "insertion within deletion",
#         "matched region within deletion",
#         "insertions within deletion and matched region",
#     ],
# )
# def test_reflect_sv_deletion_in_midsv(cssplits: str, expected: str):
#     assert cssplits_handler.reflect_sv_deletion_in_midsv(cssplits) == expected
