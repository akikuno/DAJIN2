from __future__ import annotations

from DAJIN2.core.clustering.strand_bias_handler import (
    STRAND_BIAS_LOWER_LIMIT,
    STRAND_BIAS_UPPER_LIMIT,
    count_strand,
    determine_strand_biases,
)


def test_count_strand_basic():
    labels = [1, 2, 1, 3]
    samples = [{"STRAND": "+"}, {"STRAND": "-"}, {"STRAND": "+"}, {"STRAND": "+"}]
    expected_positive = {1: 2, 3: 1}
    expected_total = {1: 2, 2: 1, 3: 1}
    positive, total = count_strand(labels, samples)
    assert positive == expected_positive
    assert total == expected_total


def test_count_strand_empty():
    labels = []
    samples = []
    expected_positive = {}
    expected_total = {}
    positive, total = count_strand(labels, samples)
    assert positive == expected_positive
    assert total == expected_total


def test_count_strand_no_positive():
    labels = [1, 1, 2]
    samples = [{"STRAND": "-"}, {"STRAND": "-"}, {"STRAND": "-"}]
    expected_positive = {}
    expected_total = {1: 2, 2: 1}
    positive, total = count_strand(labels, samples)
    assert positive == expected_positive
    assert total == expected_total


def test_determine_strand_biases():
    # Test data
    positive_counts = {1: 99, 2: 20, 3: 50}
    total_counts = {1: 100, 2: 100, 3: 100}

    expected_biases = {
        1: True,
        2: not (STRAND_BIAS_LOWER_LIMIT < 0.2 < STRAND_BIAS_UPPER_LIMIT),
        3: not (STRAND_BIAS_LOWER_LIMIT < 0.5 < STRAND_BIAS_UPPER_LIMIT),
    }

    # Call the function
    biases = determine_strand_biases(positive_counts, total_counts)

    # Assertions to check function output against expected biases
    assert biases == expected_biases


def test_edge_cases():
    # Edge cases where the ratio is exactly on the boundaries
    positive_counts = {1: 10, 2: 90}
    total_counts = {1: 100, 2: 100}

    expected_biases = {
        1: not (STRAND_BIAS_LOWER_LIMIT < 0.1 < STRAND_BIAS_UPPER_LIMIT),
        2: not (STRAND_BIAS_LOWER_LIMIT < 0.9 < STRAND_BIAS_UPPER_LIMIT),
    }

    biases = determine_strand_biases(positive_counts, total_counts)

    assert biases == expected_biases
