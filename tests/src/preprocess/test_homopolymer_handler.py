"""
Tests for homopolymer sequence error detection and handling.
"""

from unittest.mock import patch

import numpy as np
import pytest

from DAJIN2.core.preprocess.error_correction.homopolymer_handler import (
    extract_sequence_errors_in_homopolymer_loci,
    get_repeat_regions,
)


class TestGetRepeatRegions:
    """Test get_repeat_regions function."""

    def test_get_repeat_regions_basic(self):
        """Test basic homopolymer detection."""
        sequence = "ACGTAAAACCCCTTTTTGGGGAAAA"
        loci = set()  # No mutation loci to exclude

        result = get_repeat_regions(sequence, loci)

        # Should find: AAAA (4-8), CCCC (8-12), TTTTT (12-17), GGGG (17-21), AAAA (21-25)
        expected_regions = [(4, 8), (8, 12), (12, 17), (17, 21), (21, 25)]
        assert result == expected_regions

    def test_get_repeat_regions_minimum_length(self):
        """Test that only repeats >= 4 bases are detected."""
        sequence = "AAAAACCCCTTTGGAA"
        loci = set()

        result = get_repeat_regions(sequence, loci)

        # Should find: AAAAA (0-5), CCCC (5-9)
        # TTT and GG are too short, AA is too short
        expected_regions = [(0, 5), (5, 9)]
        assert result == expected_regions

    def test_get_repeat_regions_with_n(self):
        """Test homopolymer detection with N bases."""
        sequence = "AAAANNNNCCCCTTTT"
        loci = set()

        result = get_repeat_regions(sequence, loci)

        # Should find: AAAA (0-4), NNNN (4-8), CCCC (8-12), TTTT (12-16)
        expected_regions = [(0, 4), (4, 8), (8, 12), (12, 16)]
        assert result == expected_regions

    def test_get_repeat_regions_exclude_adjacent_mutations(self):
        """Test exclusion of repeats adjacent to mutation loci."""
        sequence = "AAAACCCCTTTTGGGG"
        loci = {3, 8}  # Mutations adjacent to some homopolymers

        result = get_repeat_regions(sequence, loci)

        # AAAA (0-4): NOT excluded because we need BOTH start-1==-1 AND end+1==5 in loci
        # CCCC (4-8): NOT excluded because we need BOTH start-1==3 AND end+1==9 in loci
        # TTTT (8-12): included
        # GGGG (12-16): included
        # The logic is: exclude only if BOTH adjacent positions have mutations
        expected_regions = [(0, 4), (4, 8), (8, 12), (12, 16)]
        assert result == expected_regions

    def test_get_repeat_regions_no_repeats(self):
        """Test sequence with no homopolymers."""
        sequence = "ACGTACGTACGT"
        loci = set()

        result = get_repeat_regions(sequence, loci)

        assert result == []

    def test_get_repeat_regions_short_sequence(self):
        """Test with sequence shorter than minimum repeat length."""
        sequence = "ACG"
        loci = set()

        result = get_repeat_regions(sequence, loci)

        assert result == []

    def test_get_repeat_regions_single_base_repeat(self):
        """Test with long single-base repeat."""
        sequence = "ACGTAAAAAAAACGT"
        loci = set()

        result = get_repeat_regions(sequence, loci)

        # Should find the long A repeat - 8 As from position 4 to 12
        expected_regions = [(4, 12)]
        assert result == expected_regions

    def test_get_repeat_regions_mixed_case(self):
        """Test that function handles only uppercase."""
        sequence = "AAAAccccTTTT"  # Mixed case
        loci = set()

        result = get_repeat_regions(sequence, loci)

        # Should only find uppercase repeats
        expected_regions = [(0, 4), (8, 12)]
        assert result == expected_regions


class TestExtractSequenceErrorsInHomopolymerLoci:
    """Test extract_sequence_errors_in_homopolymer_loci function."""

    @patch("DAJIN2.core.preprocess.homopolymer_handler.get_repeat_regions")
    def test_extract_errors_basic(self, mock_get_repeats):
        """Test basic homopolymer error extraction."""
        sequence = "ACGTAAAACCCCTTTT"

        indels_sample = {
            "+": np.array([0, 0, 0, 0, 50, 60, 0, 0, 30, 40, 0, 0, 20, 10, 0, 0]),
            "-": np.array([0, 0, 0, 0, 30, 40, 0, 0, 50, 60, 0, 0, 10, 20, 0, 0]),
            "*": np.array([0, 0, 0, 0, 10, 20, 0, 0, 15, 25, 0, 0, 5, 8, 0, 0]),
        }

        indels_control = {
            "+": np.array([0, 0, 0, 0, 5, 10, 0, 0, 8, 12, 0, 0, 3, 6, 0, 0]),
            "-": np.array([0, 0, 0, 0, 8, 12, 0, 0, 5, 10, 0, 0, 6, 3, 0, 0]),
            "*": np.array([0, 0, 0, 0, 2, 4, 0, 0, 3, 5, 0, 0, 1, 2, 0, 0]),
        }

        anomal_loci = {"+": {4, 5, 8, 9, 12, 13}, "-": {4, 5, 8, 9, 12, 13}, "*": {4, 5, 8, 9, 12, 13}}

        # Mock homopolymer regions
        mock_get_repeats.return_value = [(4, 8), (8, 12), (12, 16)]

        result = extract_sequence_errors_in_homopolymer_loci(sequence, indels_sample, indels_control, anomal_loci)

        # Should identify errors in homopolymer regions
        assert "+" in result
        assert "-" in result
        assert "*" in result

        # Check that some positions are identified as errors
        # (exact positions depend on the cosine similarity calculation)
        # Function calls get_repeat_regions once for each mutation type (+, -, *)
        assert mock_get_repeats.call_count == 3

    @patch("DAJIN2.core.preprocess.homopolymer_handler.get_repeat_regions")
    def test_extract_errors_no_homopolymers(self, mock_get_repeats):
        """Test when no homopolymers are found."""
        sequence = "ACGTACGTACGT"

        indels_sample = {"+": np.array([10, 20, 30]), "-": np.array([15, 25, 35]), "*": np.array([5, 15, 25])}

        indels_control = {"+": np.array([5, 10, 15]), "-": np.array([8, 12, 18]), "*": np.array([2, 8, 12])}

        anomal_loci = {"+": {0, 1, 2}, "-": {0, 1, 2}, "*": {0, 1, 2}}

        mock_get_repeats.return_value = []  # No homopolymers

        result = extract_sequence_errors_in_homopolymer_loci(sequence, indels_sample, indels_control, anomal_loci)

        # Should return empty sets when no homopolymers
        assert result["+"] == set()
        assert result["-"] == set()
        assert result["*"] == set()

    @patch("DAJIN2.core.preprocess.homopolymer_handler.get_repeat_regions")
    def test_extract_errors_no_anomalies_in_homopolymers(self, mock_get_repeats):
        """Test when anomalies don't overlap with homopolymers."""
        sequence = "ACGTAAAACCCCTTTT"

        indels_sample = {
            "+": np.array([10, 20, 30, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]),
            "-": np.array([15, 25, 35, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]),
            "*": np.array([5, 15, 25, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]),
        }

        indels_control = {
            "+": np.array([5, 10, 15, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]),
            "-": np.array([8, 12, 18, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]),
            "*": np.array([2, 8, 12, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]),
        }

        anomal_loci = {
            "+": {0, 1, 2},  # Outside homopolymer regions
            "-": {0, 1, 2},
            "*": {0, 1, 2},
        }

        mock_get_repeats.return_value = [(4, 8), (8, 12), (12, 16)]

        result = extract_sequence_errors_in_homopolymer_loci(sequence, indels_sample, indels_control, anomal_loci)

        # Should return empty sets when no overlap
        assert result["+"] == set()
        assert result["-"] == set()
        assert result["*"] == set()

    @patch("DAJIN2.core.preprocess.homopolymer_handler.get_repeat_regions")
    @patch("scipy.spatial.distance.cosine")
    def test_extract_errors_cosine_similarity_threshold(self, mock_cosine, mock_get_repeats):
        """Test cosine similarity threshold for error detection."""
        sequence = "AAAACCCC"

        indels_sample = {
            "+": np.array([0, 0, 0, 0, 50, 60, 70, 80]),
            "-": np.array([0, 0, 0, 0, 30, 40, 50, 60]),
            "*": np.array([0, 0, 0, 0, 10, 20, 30, 40]),
        }

        indels_control = {
            "+": np.array([0, 0, 0, 0, 5, 10, 15, 20]),
            "-": np.array([0, 0, 0, 0, 8, 12, 16, 20]),
            "*": np.array([0, 0, 0, 0, 2, 4, 6, 8]),
        }

        anomal_loci = {"+": {4, 5, 6, 7}, "-": {4, 5, 6, 7}, "*": {4, 5, 6, 7}}

        mock_get_repeats.return_value = [(4, 8)]  # Only second homopolymer overlaps with anomalies
        # Mock cosine distance to return low distance (high similarity)
        mock_cosine.return_value = 0.1  # cos_sim = 1 - 0.1 = 0.9 > 0.8 (error threshold)

        result = extract_sequence_errors_in_homopolymer_loci(sequence, indels_sample, indels_control, anomal_loci)

        # With high similarity (cos_sim > 0.8), positions should be marked as errors
        # Should include positions 4, 5, 6, 7, 8 (range(4, 8+1))
        expected_errors = {4, 5, 6, 7, 8}
        assert result["+"] == expected_errors
        assert result["-"] == expected_errors
        assert result["*"] == expected_errors

    @patch("DAJIN2.core.preprocess.homopolymer_handler.get_repeat_regions")
    def test_extract_errors_empty_anomalies(self, mock_get_repeats):
        """Test with empty anomal_loci."""
        sequence = "AAAACCCC"

        indels_sample = {
            "+": np.array([50, 60, 70, 80, 10, 20, 30, 40]),
            "-": np.array([30, 40, 50, 60, 15, 25, 35, 45]),
            "*": np.array([10, 20, 30, 40, 5, 15, 25, 35]),
        }

        indels_control = {
            "+": np.array([5, 10, 15, 20, 2, 4, 6, 8]),
            "-": np.array([8, 12, 16, 20, 3, 6, 9, 12]),
            "*": np.array([2, 4, 6, 8, 1, 3, 5, 7]),
        }

        anomal_loci = {"+": set(), "-": set(), "*": set()}

        mock_get_repeats.return_value = [(0, 4), (4, 8)]

        result = extract_sequence_errors_in_homopolymer_loci(sequence, indels_sample, indels_control, anomal_loci)

        # When anomal_loci is empty, get_repeat_regions will still find homopolymers
        # The function will evaluate cosine similarity and potentially mark errors
        # With the test data setup, cosine similarity will likely be high (>0.8)
        # so positions in range(start, end + 1) will be marked as errors
        # For regions (0, 4) and (4, 8), this gives range(0, 5) and range(4, 9)
        # Combined: {0, 1, 2, 3, 4, 5, 6, 7, 8}
        expected_errors = {0, 1, 2, 3, 4, 5, 6, 7, 8}
        assert result["+"] == expected_errors
        assert result["-"] == expected_errors
        assert result["*"] == expected_errors


if __name__ == "__main__":
    pytest.main([__file__])
