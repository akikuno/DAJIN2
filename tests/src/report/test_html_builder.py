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
                "</span>",
                "=G",
            ],
        ),
        # Test case 3: Consecutive Deletions
        (["=A", "-A", "-G", "=G"], ["=A", "<span class='Del_Allele'>", "-A", "-G", "</span>", "=G"]),
        # Test case 4: Single Deletion
        (["=A", "-G", "=G"], ["=A", "<span class='Del_Allele'>", "-G", "</span>", "=G"]),
        # Test case 5: Consecutive Inversions
        (["=A", "=a", "=c", "=g", "=T"], ["=A", "<span class='Inv_Allele'>", "=a", "=c", "=g", "</span>", "=T"]),
    ],
)
def test_highlight_sv_regions(misdv_sv_allele, expected_output):
    assert html_builder.highlight_sv_regions(misdv_sv_allele) == expected_output


test_data = [
    # Case 1: Insertion
    (["=A", "+A|+C|+G|=T", "=A"], ["=A", "=A", "=C", "=G", "=T", "=A"], "A<span class='Ins_Allele'>ACG</span>TA"),
    # Case 2: Insertion followed by Deletion
    (
        ["=A", "+A|+C|+G|=T", "=A"],
        ["=A", "=A", "-C", "=G", "=T", "=A"],
        "A<span class='Ins_Allele'>A<span class='Del'>C</span>G</span>TA",
    ),
    # Case 3: Insertion followed by Substitution
    (
        ["=A", "+A|+C|+G|=T", "=A"],
        ["=A", "=A", "*CG", "=G", "=T", "=A"],
        "A<span class='Ins_Allele'>A<span class='Sub'>G</span>G</span>TA",
    ),
    # Case 4: Deletion
    (["=A", "-T", "-G", "=C", "=T"], ["=A", "=C", "=T"], "A<span class='Del_Allele'>TG</span>CT"),
    # Case 5: Consecutive Deletions
    (["=A", "=T", "=C", "=T"], ["=A", "-G", "-C", "=T"], "A<span class='Del'>GC</span>T"),
    # Case 6: Complex Case with Deletions
    (
        ["=A", "-T", "-G", "=C", "=T"],
        ["=A", "-C", "=T"],
        "A<span class='Del_Allele'>TG</span><span class='Del'>C</span>T",
    ),
    # Case 7: Deletion and Substitution
    (
        ["=A", "-T", "-G", "=C", "=T"],
        ["=A", "-C", "*TG"],
        "A<span class='Del_Allele'>TG</span><span class='Del'>C</span><span class='Sub'>G</span>",
    ),
    # Case 8: Inversion
    (["=A", "=a", "=c", "=g", "=T"], ["=A", "=A", "=G", "=G", "=T"], "A<span class='Inv_Allele'>AGG</span>T"),
]


@pytest.mark.parametrize("midsv_reference, midsv_consensus, expected_output", test_data)
def test_embed_mutations_to_html(midsv_reference, midsv_consensus, expected_output):
    highlighted_sv_allele = html_builder.highlight_sv_regions(midsv_reference)
    result = html_builder.embed_mutations_to_html(highlighted_sv_allele, midsv_consensus)
    assert result == expected_output
