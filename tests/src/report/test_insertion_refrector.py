from __future__ import annotations

import pytest

from src.DAJIN2.core.report.insertion_reflector import (
    apply_insertion,
    convert_to_cstag,
    get_index_of_insertions,
    split_cstag,
)


@pytest.mark.parametrize(
    "ref_cstag, expected",
    [
        ("=AT+aatct=TGAGTTT", {2, 3, 4, 5, 6}),
        ("=AT+t=TGA+c=GTTT", {2, 6}),
        ("=ATTGAGTTT", set()),
        ("", set()),
    ],
)
def test_get_index_of_insertions(ref_cstag, expected):
    assert get_index_of_insertions(ref_cstag) == expected


@pytest.mark.parametrize(
    "cs_tag, expected",
    [
        ("=AT+t=TGA+c=GTTT", ["=A", "=T", "+t", "=T", "=G", "=A", "+c", "=G", "=T", "=T", "=T"]),
        ("=AT+t=T-ga+c=GTTT", ["=A", "=T", "+t", "=T", "-g", "-a", "+c", "=G", "=T", "=T", "=T"]),
        ("=AT+t=T*ga*ag+c=GTTT", ["=A", "=T", "+t", "=T", "*ga", "*ag", "+c", "=G", "=T", "=T", "=T"]),
        ("", []),
    ],
)
def test_split_cstag(cs_tag: str, expected: list[str]):
    assert split_cstag(cs_tag) == expected


@pytest.mark.parametrize(
    "cs_split, index_of_insertions, expected",
    [
        (
            ["=A", "=T", "+t", "=T", "=G", "=A", "+c", "=G", "=T", "=T", "=T"],
            {2, 3, 4, 5, 6},
            ["=A", "=T", "+t", "+t", "+g", "+a", "+c", "+g", "+t", "=T", "=T"],
        ),
        (
            ["=A", "=T", "+t", "=T", "-g", "=A", "+c", "=G", "=T", "=T", "=T"],
            {2, 3, 4, 5, 6},
            ["=A", "=T", "+t", "+t", "+a", "+c", "+g", "+t", "=T", "=T"],
        ),
    ],
)
def test_apply_insertion(cs_split: list[str], index_of_insertions: set[int], expected):
    assert apply_insertion(cs_split, index_of_insertions) == expected


@pytest.mark.parametrize(
    "cs_insertion, expected",
    [
        (["=A", "=T", "+t", "+t", "+a", "+c", "+g", "+t", "=T", "=T"], "=AT+ttacgt=TT"),
        (["=A", "*gt", "*gt"], "=A*gt*gt"),
        (["=A", "-g", "-g", "*gt"], "=A-gg*gt"),
    ],
)
def test_convert_to_cstag(cs_insertion: list[str], expected: str):
    assert convert_to_cstag(cs_insertion) == expected
