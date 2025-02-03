
import pytest

from DAJIN2.core.preprocess.sv_detector import get_index_and_sv_size

@pytest.mark.parametrize(
    "misdv_string, sv_type, mutation_loci, index_converter, expected",
    [
        # deletion: sv_size >= 10 のケース
        (
            "=G,-T,-T,-T,-T,-T,-T,-T,-T,-T,-T,=G",
            "deletion",
            [set()] + [{"+", "-", "*"}] * 10 + [set()],
            None,
            {1: 10},
        ),
        # deletion: sv_size < 10 のケース (戻り値は空辞書)
        (
            "=G,-T,-T,-T,-T,-T,-T,-T,-T,-T,=G",
            "deletion",
            [set()] + [{"+", "-", "*"}] * 9 + [set()],
            None,
            {},
        ),
        # inversion: sv_size >= 10 のケース
        (
            "=G,=t,=t,=t,=t,=t,=t,=t,=t,=t,=t,=G",
            "inversion",
            [set()] * 12,
            None,
            {1: 10},
        ),
        # inversion: sv_size < 10 のケース (戻り値は空辞書)
        (
            "=G,=t,=t,=t,=t,=t,=t,=t,=t,=t,=G",
            "inversion",
            [set()] * 11,
            None,
            {},
        ),
        # insertion: sv_size >= 10 のケース
        (
            "=G,+T|+T|+T|+T|+T|+T|+T|+T|+T|+T|,+T,=G",
            "insertion",
            [set(), {"+", "-"}, set(), set()],
            None,
            {1: 10},
        ),
        # insertion: sv_size < 10 のケース (戻り値は空辞書)
        (
            "=G,+T|+T|+T|+T|+T|,+T,=G",
            "insertion",
            [set(), {"+", "-"}, set(), set()],
            None,
            {},
        ),
        # index_converter の適用テスト
        (
            "=G,-T,-T,-T,-T,-T,-T,-T,-T,-T,-T,=G",
            "deletion",
            [set()] + [{"+", "-", "*"}] * 10 + [set()],
            {1: 100},  # 1 -> 100 に変換されるはず
            {100: 10},
        ),
    ],
)
def test_get_index_and_sv_size(misdv_string, sv_type, mutation_loci, index_converter, expected):
    result = get_index_and_sv_size(misdv_string, sv_type, mutation_loci, index_converter)
    assert result == expected

