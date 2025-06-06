from __future__ import annotations

import pytest

from DAJIN2.core.clustering.strand_bias_handler import (
    count_strand_by_label,
    determine_strand_biases,
)


@pytest.mark.parametrize(
    "samples, labels, expected",
    [
        # ケース1: 基本ケース
        (
            [{"STRAND": "+"}, {"STRAND": "-"}, {"STRAND": "+"}, {"STRAND": "+"}],
            [1, 2, 1, 3],
            {1: {"+": 2, "-": 0}, 2: {"+": 0, "-": 1}, 3: {"+": 1, "-": 0}},
        ),
        # ケース2: すべて+のストランド
        ([{"STRAND": "+"}, {"STRAND": "+"}], [0, 0], {0: {"+": 2, "-": 0}}),
        # ケース3: すべて-のストランド
        ([{"STRAND": "-"}, {"STRAND": "-"}], [1, 1], {1: {"+": 0, "-": 2}}),
        # ケース4: ラベルが1つずつ異なる
        ([{"STRAND": "+"}, {"STRAND": "-"}], [10, 20], {10: {"+": 1, "-": 0}, 20: {"+": 0, "-": 1}}),
        # ケース5: 空の入力
        ([], [], {}),
    ],
)
def test_count_strand_by_label(samples, labels, expected):
    assert count_strand_by_label(samples, labels) == expected


@pytest.mark.parametrize(
    "strand_counts_by_labels, strand_sample, expected",
    [
        # no bias: 比率 0.5、p-value = 1.0
        ({1: {"+": 100, "-": 100}}, {"+": 100, "-": 100}, {1: False}),
        # bias: 比率 = 0.99（偏りあり）かつ p-value 有意
        ({1: {"+": 100, "-": 1}}, {"+": 100, "-": 100}, {1: True}),
        # no bias: 統計的には有意差ありでも、比率は 0.25（偏りなしとみなす）
        ({1: {"+": 50, "-": 150}}, {"+": 150, "-": 50}, {1: False}),
    ],
)
def test_determine_strand_biases(strand_counts_by_labels, strand_sample, expected):
    assert determine_strand_biases(strand_counts_by_labels, strand_sample) == expected
