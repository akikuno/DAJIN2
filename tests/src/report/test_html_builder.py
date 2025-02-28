from __future__ import annotations

import pytest

from src.DAJIN2.core.report import html_builder


@pytest.mark.parametrize(
    "midsv_sv_allele, expected",
    [
        (["=A", "=c", "=g", "=A"], [(1, 2)]),
        (["=A", "=c", "=g", "=A", "=c", "=A"], [(1, 2), (4, 4)]),
        (["=A", "=T", "=G", "=C"], []),  # 小文字がない場合
        (["=c", "=g", "=a", "=T"], [(0, 2)]),  # すべて小文字の後に大文字
        (["=A", "=c", "=G", "=T", "=c", "=g"], [(1, 1), (4, 5)]),  # 複数の区間
        ([], []),  # 空リスト
        (["=c"], [(0, 0)]),  # 小文字1つだけ
        (["=A"], []),  # 大文字1つだけ
    ],
)
def test_get_inversion_index(midsv_sv_allele, expected):
    assert html_builder.get_inversion_index(midsv_sv_allele) == expected


@pytest.mark.parametrize(
    "midsv_sv_allele, expected",
    [
        (["=A", "+C|+C|=A", "=A"], [(1, 2)]),
        (["=A", "+C|+C|=A", "=A", "+C|+C|+C|=A", "=A"], [(1, 2), (5, 7)]),
        (["=A", "=T", "=G", "=C"], []),  # 挿入がない場合
        (["+A|+T|+G|=A", "=C"], [(0, 2)]),  # 先頭に挿入がある場合
        (["=A", "+C|+C|=A", "+T|+G|=C"], [(1, 2), (4, 5)]),  # 複数の挿入
        ([], []),  # 空リスト
        (["+A|=A"], [(0, 0)]),  # 挿入1つだけ
        (["=A"], []),  # 挿入なし
    ],
)
def test_get_insertion_index(midsv_sv_allele, expected):
    assert html_builder.get_insertion_index(midsv_sv_allele) == expected


@pytest.mark.parametrize(
    "highlight_sv_allele, expected",
    [
        (
            ["<span>=A</span>", "<span>-C", "=C", "</span>"],
            ["<span>", "=A", "</span>", "<span>", "-C", "=C", "</span>"],
        ),
        (["=A", "<b>B</b>", "C"], ["=A", "<b>", "B", "</b>", "C"]),
        (["<div><span>A</span></div>"], ["<div>", "<span>", "A", "</span>", "</div>"]),
        (["=A", "=T"], ["=A", "=T"]),  # タグなしの場合はそのまま
        (["<p>Hello <b>World</b>!</p>"], ["<p>", "Hello ", "<b>", "World", "</b>", "!", "</p>"]),
        ([], []),  # 空リスト
        (["<span>"], ["<span>"]),  # 開始タグのみ
        (["</span>"], ["</span>"]),  # 終了タグのみ
    ],
)
def test_split_html_tags(highlight_sv_allele, expected):
    assert html_builder.split_html_tags(highlight_sv_allele) == expected


###############################################################################
# append_sv_allele
###############################################################################


@pytest.mark.parametrize(
    "midsv_sv_allele, midsv_consensus, expected",
    [
        (["=A", "=c", "=g", "=A"], ["=A", "=C", "=G", "=A"], ["=A", "<span class='Inv_Allele'>=C", "=G</span>", "=A"]),
        (
            ["=A", "=c", "=g", "=A", "=c", "=A"],
            ["=A", "=C", "=G", "=A", "-C", "=A"],
            ["=A", "<span class='Inv_Allele'>=C", "=G</span>", "=A", "<span class='Inv_Allele'>-C</span>", "=A"],
        ),
    ],
)
def test_append_inversion_allele(midsv_sv_allele, midsv_consensus, expected):
    assert html_builder.append_inversion_allele(midsv_sv_allele, midsv_consensus) == expected


@pytest.mark.parametrize(
    "midsv_sv_allele, midsv_consensus, expected",
    [
        (
            ["=A", "+C|+C|=A", "=A"],
            ["=A", "=C", "=C", "=A", "=A"],
            ["=A", "<span class='Ins_Allele'>=C", "=C</span>", "=A", "=A"],
        ),
        (
            ["=A", "+C|+C|=A", "=A", "+C|+C|+C|=A", "=A"],
            ["=A", "=C", "=C", "=A", "=A", "*CG", "-C", "=C", "=A", "=A"],
            [
                "=A",
                "<span class='Ins_Allele'>=C",
                "=C</span>",
                "=A",
                "=A",
                "<span class='Ins_Allele'>*CG",
                "-C",
                "=C</span>",
                "=A",
                "=A",
            ],
        ),
    ],
)
def test_append_insertion_allele(midsv_sv_allele, midsv_consensus, expected):
    assert html_builder.append_insertion_allele(midsv_sv_allele, midsv_consensus) == expected


###############################################################################
# embed_mutations_to_html
###############################################################################


@pytest.mark.parametrize(
    "highlight_sv_allele, expected",
    [
        # Insertion
        (["=A", "+C|+C|=A", "=A"], ["<p class='p_seq'>", "A", "<span class='Ins'>CC</span>", "A", "A", "</p>"]),
        # Substitution
        (
            ["=A", "*CG", "=A", "*CG", "*CT", "=A"],
            ["<p class='p_seq'>", "A", "<span class='Sub'>G</span>", "A", "<span class='Sub'>GT</span>", "A", "</p>"],
        ),
        # Deletion
        (["=A", "-G", "-T", "=A"], ["<p class='p_seq'>", "A", "<span class='Del'>GT</span>", "A", "</p>"]),
        # Inversion
        (
            ["=A", "=c", "=g", "=A", "=c", "=A"],
            ["<p class='p_seq'>", "A", "<span class='Inv'>cg</span>", "A", "<span class='Inv'>c</span>", "A", "</p>"],
        ),
        # No mutation, just sequence
        (["=A", "=T", "=G", "=C"], ["<p class='p_seq'>", "A", "T", "G", "C", "</p>"]),
        # Empty list
        ([], ["<p class='p_seq'>", "</p>"]),
        # Combination of insertion, deletion, and substitution
        (
            ["=A", "+C|+G|=A", "-T", "*CG", "=A"],
            [
                "<p class='p_seq'>",
                "A",
                "<span class='Ins'>CG</span>",
                "A",
                "<span class='Del'>T</span>",
                "<span class='Sub'>G</span>",
                "A",
                "</p>",
            ],
        ),
    ],
)
def test_embed_mutations_to_html(highlight_sv_allele, expected):
    assert html_builder.embed_mutations_to_html(highlight_sv_allele) == expected
