import pytest

from DAJIN2.core.preprocess.structural_variants.sv_detector import get_index_and_sv_size


@pytest.mark.parametrize(
    "misdv_string, sv_type, mutation_loci, index_converter, expected",
    [
        # Deletion: sv_size >= 10
        (
            "=G,-T,-T,-T,-T,-T,-T,-T,-T,-T,-T,=G",
            "deletion",
            [set()] + [{"+", "-", "*"}] * 10 + [set()],
            None,
            {1: 10},
        ),
        # Deletion: sv_size < 10 (returns empty dict)
        (
            "=G,-T,-T,-T,-T,-T,-T,-T,-T,-T,=G",
            "deletion",
            [set()] + [{"+", "-", "*"}] * 9 + [set()],
            None,
            {},
        ),
        # Inversion: sv_size >= 10
        (
            "=G,=t,=t,=t,=t,=t,=t,=t,=t,=t,=t,=G",
            "inversion",
            [set()] * 12,
            None,
            {1: 10},
        ),
        # Inversion: sv_size < 10 (returns empty dict)
        (
            "=G,=t,=t,=t,=t,=t,=t,=t,=t,=t,=G",
            "inversion",
            [set()] * 11,
            None,
            {},
        ),
        # Insertion: sv_size >= 10
        (
            "=G,+T|+T|+T|+T|+T|+T|+T|+T|+T|+T|,+T,=G",
            "insertion",
            [set(), {"+", "-"}, set(), set()],
            None,
            {1: 10},
        ),
        # Insertion: sv_size < 10 (returns empty dict)
        (
            "=G,+T|+T|+T|+T|+T|,+T,=G",
            "insertion",
            [set(), {"+", "-"}, set(), set()],
            None,
            {},
        ),
        # Test index_converter application
        (
            "=G,-T,-T,-T,-T,-T,-T,-T,-T,-T,-T,=G",
            "deletion",
            [set()] + [{"+", "-", "*"}] * 10 + [set()],
            {1: 100},  # 1 should be mapped to 100
            {100: 10},
        ),
    ],
)
def test_get_index_and_sv_size(misdv_string, sv_type, mutation_loci, index_converter, expected):
    result = get_index_and_sv_size(misdv_string, sv_type, mutation_loci, index_converter)
    assert result == expected
