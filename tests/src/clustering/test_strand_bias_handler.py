from __future__ import annotations

import pytest

from DAJIN2.core.clustering.strand_bias_handler import (
    count_strand_by_label,
    determine_strand_biases,
)


@pytest.mark.parametrize(
    "samples, labels, expected",
    [
        # Case 1: Basic case
        (
            [{"STRAND": "+"}, {"STRAND": "-"}, {"STRAND": "+"}, {"STRAND": "+"}],
            [1, 2, 1, 3],
            {1: {"+": 2, "-": 0}, 2: {"+": 0, "-": 1}, 3: {"+": 1, "-": 0}},
        ),
        # Case 2: All plus strand
        ([{"STRAND": "+"}, {"STRAND": "+"}], [0, 0], {0: {"+": 2, "-": 0}}),
        # Case 3: All minus strand
        ([{"STRAND": "-"}, {"STRAND": "-"}], [1, 1], {1: {"+": 0, "-": 2}}),
        # Case 4: Each label is different
        ([{"STRAND": "+"}, {"STRAND": "-"}], [10, 20], {10: {"+": 1, "-": 0}, 20: {"+": 0, "-": 1}}),
        # Case 5: Empty input
        ([], [], {}),
    ],
)
def test_count_strand_by_label(samples, labels, expected):
    assert count_strand_by_label(samples, labels) == expected


@pytest.mark.parametrize(
    "strand_counts_by_labels, strand_sample, expected",
    [
        # No bias: ratio 0.5, p-value = 1.0
        ({1: {"+": 100, "-": 100}}, {"+": 100, "-": 100}, {1: False}),
        # Bias: ratio = 0.99 (skewed) and p-value is significant
        ({1: {"+": 100, "-": 1}}, {"+": 100, "-": 100}, {1: True}),
        # No bias: statistically significant, but ratio 0.25 is treated as not skewed
        ({1: {"+": 50, "-": 150}}, {"+": 150, "-": 50}, {1: False}),
    ],
)
def test_determine_strand_biases(strand_counts_by_labels, strand_sample, expected):
    assert determine_strand_biases(strand_counts_by_labels, strand_sample) == expected
