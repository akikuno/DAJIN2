import pytest

from DAJIN2.src.DAJIN2.core.consensus import sv_annotator


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
