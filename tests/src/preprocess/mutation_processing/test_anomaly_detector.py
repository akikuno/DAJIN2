"""
Tests for anomaly detection utilities for mutation extraction.
"""

from unittest.mock import Mock, patch

import numpy as np
import pytest

from DAJIN2.core.preprocess.mutation_processing.anomaly_detector import (
    cosine_distance,
    detect_anomalies,
    extract_anomal_loci,
    is_dissimilar_loci,
)


class TestCosineDistance:
    """Test cosine_distance function."""

    def test_cosine_distance_identical_vectors(self):
        """Test cosine distance between identical vectors."""
        x = [1.0, 2.0, 3.0]
        y = [1.0, 2.0, 3.0]

        result = cosine_distance(x, y)

        assert result == pytest.approx(0.0, abs=1e-6)

    def test_cosine_distance_orthogonal_vectors(self):
        """Test cosine distance between orthogonal vectors."""
        x = [1.0, 0.0]
        y = [0.0, 1.0]

        result = cosine_distance(x, y)

        # Due to 1e-6 addition for numerical stability, result is slightly less than 1.0
        assert result == pytest.approx(1.0, abs=2e-6)

    def test_cosine_distance_opposite_vectors(self):
        """Test cosine distance between opposite vectors."""
        x = [1.0, 1.0]
        y = [-1.0, -1.0]

        result = cosine_distance(x, y)

        assert result == pytest.approx(2.0, abs=1e-6)

    def test_cosine_distance_zero_vectors(self):
        """Test cosine distance with zero vectors (handles division by zero)."""
        x = [0.0, 0.0]
        y = [0.0, 0.0]

        # Should not raise an error due to 1e-6 addition
        result = cosine_distance(x, y)

        assert isinstance(result, float)

    def test_cosine_distance_mixed_values(self):
        """Test cosine distance with realistic values."""
        x = [10.5, 20.2, 5.8]
        y = [12.1, 18.7, 6.2]

        result = cosine_distance(x, y)

        assert 0.0 <= result <= 2.0  # Valid range for cosine distance


class TestIsDissimilarLoci:
    """Test is_dissimilar_loci function."""

    def test_is_dissimilar_loci_large_difference_consensus(self):
        """Test with large difference in consensus mode."""
        values_sample = np.array([100.0, 50.0, 80.0, 25.0, 40.0, 30.0, 20.0, 10.0, 5.0, 2.0])
        values_control = np.array([10.0, 20.0, 30.0, 5.0, 15.0, 12.0, 8.0, 4.0, 2.0, 1.0])
        index = 0

        result = is_dissimilar_loci(values_sample, values_control, index, is_consensus=True)

        # Large difference (100 - 10 = 90 > 20) and sample > 50
        assert result is True

    def test_is_dissimilar_loci_large_difference_consensus_low_sample(self):
        """Test with large difference but low sample value in consensus mode."""
        values_sample = np.array([40.0, 30.0, 25.0, 20.0, 15.0, 10.0, 8.0, 5.0, 3.0, 1.0])
        values_control = np.array([10.0, 15.0, 12.0, 8.0, 6.0, 4.0, 3.0, 2.0, 1.0, 0.5])
        index = 0

        result = is_dissimilar_loci(values_sample, values_control, index, is_consensus=True)

        # Large difference (40 - 10 = 30 > 20) but sample <= 50
        assert result is False

    def test_is_dissimilar_loci_large_difference_non_consensus(self):
        """Test with large difference in non-consensus mode."""
        values_sample = np.array([40.0, 30.0, 25.0, 20.0, 15.0, 10.0, 8.0, 5.0, 3.0, 1.0])
        values_control = np.array([10.0, 15.0, 12.0, 8.0, 6.0, 4.0, 3.0, 2.0, 1.0, 0.5])
        index = 0

        result = is_dissimilar_loci(values_sample, values_control, index, is_consensus=False)

        # Large difference regardless of sample value in non-consensus mode
        assert result is True

    def test_is_dissimilar_loci_small_difference(self):
        """Test with small difference."""
        values_sample = np.array([15.0, 12.0, 10.0, 8.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0])
        values_control = np.array([10.0, 8.0, 6.0, 5.0, 4.0, 3.0, 2.5, 2.0, 1.5, 1.0])
        index = 0

        # Mock cosine_distance to return predictable values
        with patch("DAJIN2.core.preprocess.mutation_processing.anomaly_detector.cosine_distance") as mock_cosine:
            mock_cosine.side_effect = [0.1, 0.05]  # distance, distance_slice

            result = is_dissimilar_loci(values_sample, values_control, index)

            # distance (0.1) > 0.05 and distance / (distance + distance_slice) = 0.1 / 0.15 ≈ 0.67 < 0.9
            assert result is False

    def test_is_dissimilar_loci_high_cosine_distance(self):
        """Test with high cosine distance indicating dissimilarity."""
        values_sample = np.array([15.0, 12.0, 10.0, 8.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0])
        values_control = np.array([10.0, 8.0, 6.0, 5.0, 4.0, 3.0, 2.5, 2.0, 1.5, 1.0])
        index = 0

        with patch("DAJIN2.core.preprocess.mutation_processing.anomaly_detector.cosine_distance") as mock_cosine:
            mock_cosine.side_effect = [0.1, 0.01]  # distance, distance_slice

            result = is_dissimilar_loci(values_sample, values_control, index)

            # distance (0.1) > 0.05 and distance / (distance + distance_slice) = 0.1 / 0.11 ≈ 0.91 > 0.9
            assert result is True

    def test_is_dissimilar_loci_insufficient_length(self):
        """Test with arrays shorter than required window."""
        values_sample = np.array([15.0, 12.0])
        values_control = np.array([10.0, 8.0])
        index = 0

        with patch("DAJIN2.core.preprocess.mutation_processing.anomaly_detector.cosine_distance") as mock_cosine:
            mock_cosine.side_effect = [0.1, 0.01]

            result = is_dissimilar_loci(values_sample, values_control, index)

            # Should handle short arrays gracefully
            assert isinstance(result, bool)


class TestDetectAnomalies:
    """Test detect_anomalies function."""

    @patch("DAJIN2.core.preprocess.mutation_processing.anomaly_detector.is_dissimilar_loci")
    def test_detect_anomalies_basic(self, mock_is_dissimilar):
        """Test basic anomaly detection."""
        values_sample = np.array([10.0, 20.0, 30.0, 40.0, 50.0])
        values_control = np.array([5.0, 10.0, 15.0, 20.0, 25.0])
        threshold = 5.0

        # Mock MLPClassifier to return predictable results
        with patch("DAJIN2.core.preprocess.mutation_processing.anomaly_detector.MLPClassifier") as mock_mlp:
            mock_classifier = Mock()
            mock_classifier.predict.return_value = [0, 1, 1, 0, 1]  # Predict anomalies at indices 1, 2, 4
            mock_mlp.return_value = mock_classifier

            # Mock is_dissimilar_loci to return True for indices 1 and 4
            mock_is_dissimilar.side_effect = lambda vs, vc, i, ic: i in [1, 4]

            result = detect_anomalies(values_sample, values_control, threshold)

            assert result == {1, 4}

    @patch("DAJIN2.core.preprocess.mutation_processing.anomaly_detector.is_dissimilar_loci")
    def test_detect_anomalies_no_anomalies(self, mock_is_dissimilar):
        """Test when no anomalies are detected."""
        values_sample = np.array([10.0, 20.0, 30.0])
        values_control = np.array([8.0, 18.0, 28.0])
        threshold = 5.0

        with patch("DAJIN2.core.preprocess.mutation_processing.anomaly_detector.MLPClassifier") as mock_mlp:
            mock_classifier = Mock()
            mock_classifier.predict.return_value = [0, 0, 0]  # No anomalies predicted
            mock_mlp.return_value = mock_classifier

            result = detect_anomalies(values_sample, values_control, threshold)

            assert result == set()

    @patch("DAJIN2.core.preprocess.mutation_processing.anomaly_detector.is_dissimilar_loci")
    def test_detect_anomalies_consensus_mode(self, mock_is_dissimilar):
        """Test anomaly detection in consensus mode."""
        values_sample = np.array([10.0, 20.0])
        values_control = np.array([5.0, 10.0])
        threshold = 3.0

        with patch("DAJIN2.core.preprocess.mutation_processing.anomaly_detector.MLPClassifier") as mock_mlp:
            mock_classifier = Mock()
            mock_classifier.predict.return_value = [1, 1]
            mock_mlp.return_value = mock_classifier

            mock_is_dissimilar.side_effect = lambda vs, vc, i, ic: ic  # Only return True in consensus mode

            result = detect_anomalies(values_sample, values_control, threshold, is_consensus=True)

            # In consensus mode with mock returning True for all indices
            assert result == {0, 1}

            # Verify is_dissimilar_loci was called with is_consensus=True
            for call in mock_is_dissimilar.call_args_list:
                # is_consensus is the 4th argument (index 3)
                assert call[0][3] is True

    def test_detect_anomalies_deterministic(self):
        """Test that detect_anomalies produces deterministic results."""
        values_sample = np.array([10.0, 20.0, 30.0])
        values_control = np.array([5.0, 10.0, 15.0])
        threshold = 5.0

        result1 = detect_anomalies(values_sample, values_control, threshold)
        result2 = detect_anomalies(values_sample, values_control, threshold)

        # Should be deterministic due to fixed random seed
        assert result1 == result2


class TestExtractAnomalLoci:
    """Test extract_anomal_loci function."""

    @patch("DAJIN2.core.preprocess.mutation_processing.anomaly_detector.detect_anomalies")
    def test_extract_anomal_loci_basic(self, mock_detect):
        """Test basic anomalous loci extraction."""
        indels_normalized_sample = {
            "+": np.array([10.0, 20.0, 30.0]),
            "-": np.array([5.0, 15.0, 25.0]),
            "*": np.array([2.0, 8.0, 12.0]),
        }
        indels_normalized_control = {
            "+": np.array([5.0, 10.0, 15.0]),
            "-": np.array([3.0, 8.0, 12.0]),
            "*": np.array([1.0, 4.0, 6.0]),
        }
        thresholds = {"+": 5.0, "-": 3.0, "*": 2.0}

        # Mock detect_anomalies to return different results for each mutation type
        def side_effect(vs, vc, thresh, is_consensus=False):
            # Return different results based on threshold
            if thresh == 5.0:  # insertions "+"
                return {0, 2}
            elif thresh == 3.0:  # deletions "-"
                return {1}
            elif thresh == 2.0:  # substitutions "*"
                return {0, 1, 2}
            return set()

        mock_detect.side_effect = side_effect

        result = extract_anomal_loci(indels_normalized_sample, indels_normalized_control, thresholds)

        assert result["+"] == {0, 2}
        assert result["-"] == {1}
        assert result["*"] == {0, 1, 2}

        # Verify detect_anomalies was called for each mutation type
        assert mock_detect.call_count == 3

    @patch("DAJIN2.core.preprocess.mutation_processing.anomaly_detector.detect_anomalies")
    def test_extract_anomal_loci_consensus_mode(self, mock_detect):
        """Test extraction in consensus mode."""
        indels_normalized_sample = {"+": np.array([10.0]), "-": np.array([5.0]), "*": np.array([2.0])}
        indels_normalized_control = {"+": np.array([5.0]), "-": np.array([3.0]), "*": np.array([1.0])}
        thresholds = {"+": 5.0, "-": 3.0, "*": 2.0}

        mock_detect.return_value = set()

        extract_anomal_loci(indels_normalized_sample, indels_normalized_control, thresholds, is_consensus=True)

        # Verify all calls were made with is_consensus=True
        for call in mock_detect.call_args_list:
            # is_consensus is passed as 4th positional argument (index 3)
            assert call[0][3] is True

    @patch("DAJIN2.core.preprocess.mutation_processing.anomaly_detector.detect_anomalies")
    def test_extract_anomal_loci_empty_results(self, mock_detect):
        """Test extraction with no anomalies found."""
        indels_normalized_sample = {"+": np.array([10.0, 20.0]), "-": np.array([5.0, 15.0]), "*": np.array([2.0, 8.0])}
        indels_normalized_control = {"+": np.array([8.0, 18.0]), "-": np.array([4.0, 13.0]), "*": np.array([1.5, 7.0])}
        thresholds = {"+": 5.0, "-": 3.0, "*": 2.0}

        mock_detect.return_value = set()

        result = extract_anomal_loci(indels_normalized_sample, indels_normalized_control, thresholds)

        assert result["+"] == set()
        assert result["-"] == set()
        assert result["*"] == set()


if __name__ == "__main__":
    pytest.main([__file__])
