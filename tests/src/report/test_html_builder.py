from __future__ import annotations

import pytest

from src.DAJIN2.core.report import html_builder


# Test function
@pytest.mark.parametrize(
    "misdv_sv_allele, expected_output",
    [
        # Test case 1: Insertion
        (["=A", "+A|+A|=C", "=G"], ["=A", "<span class='Ins_Allele'>", "=A", "=A", "</span>", "=C", "=G"]),
        # Test case 2: Insertion followed by Deletion
        (
            ["=A", "+A|+A|-C", "=G"],
            [
                "=A",
                "<span class='Ins_Allele'>",
                "=A",
                "=A",
                "</span>",
                "<span class='Del_Allele'>",
                "-C",
                "</span><!-- END_OF_DEL_ALLELE -->",
                "=G",
            ],
        ),
        # Test case 3: Consecutive Deletions
        (["=A", "-A", "-G", "=G"], ["=A", "<span class='Del_Allele'>", "-A", "-G", "</span><!-- END_OF_DEL_ALLELE -->", "=G"]),
        # Test case 4: Single Deletion
        (["=A", "-G", "=G"], ["=A", "<span class='Del_Allele'>", "-G", "</span><!-- END_OF_DEL_ALLELE -->", "=G"]),
        # Test case 5: Consecutive Inversions
        (["=A", "=a", "=c", "=g", "=T"], ["=A", "<span class='Inv_Allele'>", "=a", "=c", "=g", "</span>", "=T"]),
    ],
)
def test_highlight_sv_regions(misdv_sv_allele, expected_output):
    assert html_builder.highlight_sv_regions(misdv_sv_allele) == expected_output


test_data = [
    # Case 1: Insertion
    (
        ["=A", "+A|+C|+G|=T", "=A"],
        ["=A", "=A", "=C", "=G", "=T", "=A"],
        "<p class='p_seq'>A<span class='Ins_Allele'>ACG</span>TA</p>",
    ),
    # Case 2: Insertion followed by Deletion
    (
        ["=A", "+A|+C|+G|=T", "=A"],
        ["=A", "=A", "-C", "=G", "=T", "=A"],
        "<p class='p_seq'>A<span class='Ins_Allele'>A<span class='Del'>C</span>G</span>TA</p>",
    ),
    # Case 3: Insertion followed by Substitution
    (
        ["=A", "+A|+C|+G|=T", "=A"],
        ["=A", "=A", "*CG", "=G", "=T", "=A"],
        "<p class='p_seq'>A<span class='Ins_Allele'>A<span class='Sub'>G</span>G</span>TA</p>",
    ),
    # Case 4: Deletion
    (["=A", "-T", "-G", "=C", "=T"], ["=A", "=C", "=T"], "<p class='p_seq'>A<span class='Del_Allele'>TG</span><!-- END_OF_DEL_ALLELE -->CT</p>"),
    # Case 5: Consecutive Deletions
    (["=A", "=T", "=C", "=T"], ["=A", "-G", "-C", "=T"], "<p class='p_seq'>A<span class='Del'>GC</span>T</p>"),
    # Case 6: Complex Case with Deletions
    (
        ["=A", "-T", "-G", "=C", "=T"],
        ["=A", "-C", "=T"],
        "<p class='p_seq'>A<span class='Del_Allele'>TG</span><!-- END_OF_DEL_ALLELE --><span class='Del'>C</span>T</p>",
    ),
    # Case 7: Deletion and Substitution
    (
        ["=A", "-T", "-G", "=C", "=T"],
        ["=A", "-C", "*TG"],
        "<p class='p_seq'>A<span class='Del_Allele'>TG</span><!-- END_OF_DEL_ALLELE --><span class='Del'>C</span><span class='Sub'>G</span></p>",
    ),
    # Case 8: Inversion
    (
        ["=A", "=a", "=c", "=g", "=T"],
        ["=A", "=A", "=G", "=G", "=T"],
        "<p class='p_seq'>A<span class='Inv_Allele'>AGG</span>T</p>",
    ),
]


@pytest.mark.parametrize("midsv_reference, midsv_consensus, expected_output", test_data)
def test_embed_mutations_to_html(midsv_reference, midsv_consensus, expected_output):
    highlighted_sv_allele = html_builder.highlight_sv_regions(midsv_reference)
    result = html_builder.embed_mutations_to_html(highlighted_sv_allele, midsv_consensus)
    assert "".join(result) == expected_output


###############################################################################
# Handle insertions
###############################################################################

@pytest.mark.parametrize(
    "html_parts, expected",
    [
        # 提示されたテストケース 1
        (
            [
                "<span class='Del_Allele'>", "A", "</span><!-- END_OF_DEL_ALLELE -->",
                "<span>N", "<span>N", "<span>N", "<span>N", "<span>N",
                "<span class='Del_Allele'>", "A", "</span><!-- END_OF_DEL_ALLELE -->"
            ],
            (
                [['<span>N', '<span>N', '<span>N', '<span>N', '<span>N']],
                [(3, 8)]
            )
        ),
        # 提示されたテストケース 2
        (
            [
                "<span class='Del_Allele'>", "A", "</span><!-- END_OF_DEL_ALLELE -->",
                "<span>N", "<span>N", "<span>N", "<span>N", "<span>N",
                "<span class='Del_Allele'>", "A", "</span><!-- END_OF_DEL_ALLELE -->",
                "<span>T", "<span>T", "<span>T", "<span>T", "<span>T",
                "<span class='Del_Allele'>", "C", "</span><!-- END_OF_DEL_ALLELE -->"
            ],
            (
                [
                    ['<span>N', '<span>N', '<span>N', '<span>N', '<span>N'],
                    ['<span>T', '<span>T', '<span>T', '<span>T', '<span>T']
                ],
                [(3, 8), (11, 16)]
            )
        ),
        # 空リスト
        ([], ([], [])),
        # `<!-- END_OF_DEL_ALLELE -->` がない場合
        (
            ["<span class='Del_Allele'>", "A", "</span>"],
            ([], [])
        ),
        # `Del_Allele` がない場合
        (
            ["A", "B", "C", "<!-- END_OF_DEL_ALLELE -->", "D", "E"],
            ([], [])
        ),
        # `Del_Allele` が `END_OF_DEL_ALLELE` の後にこない場合
        (
            [
                "<span class='Del_Allele'>", "A", "</span><!-- END_OF_DEL_ALLELE -->",
                "B", "C", "D"
            ],
            ([], [])
        ),
    ],
)
def test_get_flanked_tags_by_deletions(html_parts, expected):
    assert html_builder.get_flanked_tags_by_deletions(html_parts) == expected


@pytest.mark.parametrize(
    "flanked_tag, n, expected",
    [
        # 10 文字以上の連続マッチがある場合（デフォルト n=10）
        (["A", "A", "A", "A", "A", "A", "A", "A", "A", "A"], 10, True),
        (["G", "G", "G", "G", "G", "G", "G", "G", "G", "G", "G"], 10, True),

        # 9 文字以下の連続マッチ（n=10 で False）
        (["T", "T", "T", "T", "T", "T", "T", "T", "T"], 10, False),
        (["C", "C", "C", "C", "C", "C", "C", "C"], 10, False),

        # HTMLタグが途中に入るとカウントがリセットされる
        (["A", "A", "<span>", "A", "A", "A", "A", "A", "A", "A", "A"], 10, False),
        (["G", "G", "G", "<b>", "G", "G", "G", "G", "G", "G", "G", "G"], 10, False),

        # 連続マッチがあるが、閾値 n を変更して判定
        (["A", "A", "A", "A", "A"], 5, True),
        (["A", "A", "A", "A"], 5, False),

        # すべて HTML タグで構成されている場合
        (["<span>", "<b>", "<div>"], 10, False),

        # 空リスト
        ([], 10, False),
    ],
)
def test_has_consecutive_matches(flanked_tag, n, expected):
    assert html_builder.has_consecutive_matches(flanked_tag, n) == expected


@pytest.mark.parametrize(
    "html_parts, sequences_within_deletion, expected",
    [
        # 0.提示されたテストケース
        (
            [
                "<span class='Del_Allele'>", "A", "</span><!-- END_OF_DEL_ALLELE -->",
                "<span>N", "<span>N", "<span>N", "<span>N", "<span>N",
                "<span class='Del_Allele'>", "A", "</span><!-- END_OF_DEL_ALLELE -->"
            ],
            [["T", "T", "T", "T", "T"]],
            [
                "<span class='Del_Allele'>", "A", 'T', 'T', 'T', 'T', 'T', 'A',
                "</span><!-- END_OF_DEL_ALLELE -->", "<span class='Ins'>NNNNN</span>"
            ]
        ),
        # 2.2つの異なるDeletion領域がある場合
        (
            [
                "<span class='Del_Allele'>", "A", "A", "</span><!-- END_OF_DEL_ALLELE -->",
                "<span>N", "<span>N",
                "<span class='Del_Allele'>", "A", "A", "</span><!-- END_OF_DEL_ALLELE -->",
                "<span>N", "<span>N",
                "<span class='Del_Allele'>", "A", "A", "</span><!-- END_OF_DEL_ALLELE -->"
            ],
            [["T","T"], ["T","T"]],
            ["<span class='Del_Allele'>", 'A', 'A', 'T', 'T', 'A', 'A', 'T', 'T', 'A', 'A', '</span><!-- END_OF_DEL_ALLELE -->', "<span class='Ins'>NNNN</span>"]
        ),
        # 3.2つの異なるDeletion領域があり、最初だけ挿入が発生する場合
        (
            [
                "<span class='Del_Allele'>", "A", "A", "</span><!-- END_OF_DEL_ALLELE -->",
                "<span>N", "<span>N",
                "<span class='Del_Allele'>", "A", "A", "</span><!-- END_OF_DEL_ALLELE -->",
                "G","G","G","G","G","G","G","G","G","G","G","G","G","G","G",
                "<span class='Del_Allele'>", "A", "A", "</span><!-- END_OF_DEL_ALLELE -->"
            ],
            [["T","T"], ["G","G","G","G","G","G","G","G","G","G","G","G","G","G","G"]],
            ["<span class='Del_Allele'>", 'A', 'A', 'T', 'T', 'A', 'A', '</span><!-- END_OF_DEL_ALLELE -->', "<span class='Ins'>NN</span>", 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', "<span class='Del_Allele'>", 'A', 'A', '</span><!-- END_OF_DEL_ALLELE -->']
        ),

        # 3.Deletion領域がなく挿入も発生しない場合
        (
            ["A", "T", "G", "C"],
            [],
            ["A", "T", "G", "C"]
        ),
    ],
)
def test_handle_insertions(html_parts, sequences_within_deletion, expected):
    flanked_tags, flanked_indices = html_builder.get_flanked_tags_by_deletions(html_parts)
    result = html_builder.handle_insertions(html_parts[:], flanked_tags, flanked_indices, sequences_within_deletion)
    assert result == expected
