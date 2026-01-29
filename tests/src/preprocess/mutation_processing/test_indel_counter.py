"""
Tests for indel counting and normalization utilities.
"""

from pathlib import Path
from unittest.mock import patch

import numpy as np
import pytest

from DAJIN2.core.preprocess.mutation_processing.indel_counter import (
    count_indels,
    minimize_mutation_counts,
    normalize_indels,
    summarize_indels,
)


class TestCountIndels:
    """Test count_indels function."""

    def test_count_indels_basic(self):
        """Test basic indel counting functionality."""
        sequence = "ACGTACGT"
        midsv_data = [
            {"MIDSV": "=A,=C,=G,=T,=A,=C,=G,=T"},
            {"MIDSV": "=A,=C,+G,=T,=A,=C,=G,=T"},
            {"MIDSV": "=A,=C,=G,-T,=A,=C,=G,=T"},
            {"MIDSV": "=A,=C,=G,*T,=A,=C,=G,=T"},
        ]

        result = count_indels(iter(midsv_data), sequence)

        assert len(result) == 4
        assert "=" in result
        assert "+" in result
        assert "-" in result
        assert "*" in result
        assert len(result["="]) == 8
        assert result["="][0] == 4  # All samples have match at position 0
        assert result["+"][2] == 1  # One insertion at position 2
        assert result["-"][3] == 1  # One deletion at position 3
        assert result["*"][3] == 1  # One substitution at position 3

    def test_count_indels_with_n_and_lowercase(self):
        """Test that N and lowercase characters are skipped."""
        sequence = "ACGT"
        midsv_data = [
            {"MIDSV": "=A,=N,=G,=T"},
            {"MIDSV": "=A,=c,=G,=T"},  # lowercase
        ]

        result = count_indels(iter(midsv_data), sequence)

        assert result["="][1] == 0  # N and lowercase are skipped
        assert result["="][0] == 2  # Both samples counted at position 0

    def test_count_indels_empty_input(self):
        """Test with empty input."""
        sequence = "ACGT"
        midsv_data = []

        result = count_indels(iter(midsv_data), sequence)

        assert all(count == 0 for counts in result.values() for count in counts)


class TestNormalizeIndels:
    """Test normalize_indels function."""

    def test_normalize_indels_basic(self):
        """Test basic normalization functionality."""
        count = {"=": [10, 8, 6, 4], "+": [0, 2, 0, 1], "-": [1, 0, 2, 0], "*": [0, 1, 1, 2]}

        result = normalize_indels(count)

        # Check that percentages are calculated correctly
        np.testing.assert_array_almost_equal(result["+"], [0.0, 20.0, 0.0, 20.0])
        np.testing.assert_array_almost_equal(result["-"], [9.09090909, 0.0, 25.0, 0.0])
        np.testing.assert_array_almost_equal(result["*"], [0.0, 11.11111111, 14.28571429, 33.33333333])

    def test_normalize_indels_zero_coverage(self):
        """Test normalization with zero coverage positions."""
        count = {"=": [0, 5], "+": [0, 1], "-": [1, 0], "*": [0, 0]}

        result = normalize_indels(count)

        # When denominator is 0, result should be 0
        # For first position: matches=0, mutations=0, so 0/(0+0) should be 0
        # For "+": 0/(0+0) = 0
        # For "-": 1/(1+0) = 100
        # For "*": 0/(0+0) = 0
        assert result["+"][0] == 0.0
        assert result["-"][0] == 100.0  # 1/(1+0)*100 = 100
        assert result["*"][0] == 0.0

    def test_normalize_indels_all_mutations(self):
        """Test normalization when all reads have mutations."""
        count = {"=": [0, 0, 0], "+": [5, 0, 0], "-": [0, 3, 0], "*": [0, 0, 2]}

        result = normalize_indels(count)

        assert result["+"][0] == 100.0
        assert result["-"][1] == 100.0
        assert result["*"][2] == 100.0


class TestMinimizeMutationCounts:
    """Test minimize_mutation_counts function."""

    def test_minimize_mutation_counts_basic(self):
        """Test basic minimization functionality."""
        indels_control = {
            "+": np.array([10.0, 20.0, 5.0]),
            "-": np.array([15.0, 5.0, 25.0]),
            "*": np.array([8.0, 12.0, 3.0]),
        }
        indels_sample = {
            "+": np.array([8.0, 25.0, 10.0]),
            "-": np.array([20.0, 3.0, 15.0]),
            "*": np.array([10.0, 8.0, 5.0]),
        }

        result = minimize_mutation_counts(indels_control, indels_sample)

        # Should take minimum between control and sample
        np.testing.assert_array_equal(result["+"], [8.0, 20.0, 5.0])
        np.testing.assert_array_equal(result["-"], [15.0, 3.0, 15.0])
        np.testing.assert_array_equal(result["*"], [8.0, 8.0, 3.0])

    def test_minimize_mutation_counts_identical(self):
        """Test when control and sample are identical."""
        indels_control = {"+": np.array([10.0, 20.0]), "-": np.array([15.0, 5.0]), "*": np.array([8.0, 12.0])}
        indels_sample = indels_control.copy()

        result = minimize_mutation_counts(indels_control, indels_sample)

        # Should return identical arrays
        np.testing.assert_array_equal(result["+"], indels_control["+"])
        np.testing.assert_array_equal(result["-"], indels_control["-"])
        np.testing.assert_array_equal(result["*"], indels_control["*"])


class TestSummarizeIndels:
    """Test summarize_indels function."""

    @patch("DAJIN2.core.preprocess.mutation_processing.indel_counter.fileio.read_jsonl")
    def test_summarize_indels_integration(self, mock_read_jsonl):
        """Test integration of count_indels and normalize_indels."""
        mock_read_jsonl.return_value = [
            {"MIDSV": "=A,=C,=G,=T"},
            {"MIDSV": "=A,+C,=G,=T"},
        ]

        path_midsv = Path("dummy.jsonl")
        sequence = "ACGT"

        indels_count, indels_normalized = summarize_indels(path_midsv, sequence)

        mock_read_jsonl.assert_called_once_with(path_midsv)

        # Check that we get both count and normalized results
        assert "=" in indels_count
        assert "+" in indels_count
        assert "=" in indels_normalized
        assert "+" in indels_normalized

        # Verify counts
        assert indels_count["="][0] == 2  # Both samples match at position 0
        assert indels_count["+"][1] == 1  # One insertion at position 1

        # Verify normalization
        assert indels_normalized["+"][1] == 50.0  # 1/2 * 100 = 50%

    @patch("DAJIN2.core.preprocess.mutation_processing.indel_counter.fileio.read_jsonl")
    def test_summarize_indels_empty_file(self, mock_read_jsonl):
        """Test summarize_indels with empty input file."""
        mock_read_jsonl.return_value = []

        path_midsv = Path("empty.jsonl")
        sequence = "ACGT"

        indels_count, indels_normalized = summarize_indels(path_midsv, sequence)

        # Should handle empty input gracefully
        assert all(count == 0 for counts in indels_count.values() for count in counts)
        assert all(count == 0.0 for counts in indels_normalized.values() for count in counts)


if __name__ == "__main__":
    pytest.main([__file__])
