"""
Tests for indel merging and mutation loci handling utilities.
"""

import numpy as np
import pytest

from DAJIN2.core.preprocess.mutation_processing.indel_merger import (
    add_knockin_loci,
    count_elements_within_range,
    discard_errors,
    merge_index_of_consecutive_indel,
    merge_loci,
    split_kmer,
    transpose_mutation_loci,
)


class TestCountElementsWithinRange:
    """Test count_elements_within_range function."""

    def test_count_elements_basic(self):
        """Test basic counting functionality."""
        arr = [1, 3, 5, 7, 9, 11, 13, 15]

        # Count elements between 5 and 11 (inclusive)
        result = count_elements_within_range(arr, 5, 11)

        assert result == 4  # [5, 7, 9, 11]

    def test_count_elements_empty_range(self):
        """Test with range containing no elements."""
        arr = [1, 3, 5, 7, 9]

        result = count_elements_within_range(arr, 10, 15)

        assert result == 0

    def test_count_elements_single_element(self):
        """Test with range containing single element."""
        arr = [1, 3, 5, 7, 9]

        result = count_elements_within_range(arr, 5, 5)

        assert result == 1

    def test_count_elements_full_range(self):
        """Test with range containing all elements."""
        arr = [2, 4, 6, 8]

        result = count_elements_within_range(arr, 0, 10)

        assert result == 4

    def test_count_elements_empty_array(self):
        """Test with empty array."""
        arr = []

        result = count_elements_within_range(arr, 1, 5)

        assert result == 0

    def test_count_elements_duplicates(self):
        """Test with duplicate values."""
        arr = [1, 2, 2, 3, 3, 3, 4]

        result = count_elements_within_range(arr, 2, 3)

        assert result == 5  # [2, 2, 3, 3, 3]


class TestMergeIndexOfConsecutiveIndel:
    """Test merge_index_of_consecutive_indel function."""

    def test_merge_consecutive_basic(self):
        """Test basic consecutive indel merging."""
        mutation_loci = {"+": {5, 10, 15}, "-": {20, 25, 30}, "*": {1, 2, 3}}

        result = merge_index_of_consecutive_indel(mutation_loci)

        # Point mutations should remain unchanged
        assert result["*"] == {1, 2, 3}

        # Indels should be checked for merging
        assert "+" in result
        assert "-" in result

    def test_merge_consecutive_close_indels(self):
        """Test merging of close indels (within 10 bases)."""
        mutation_loci = {
            "+": {5, 12},  # 12 - 5 = 7 < 10, should merge
            "-": {20, 35},  # 35 - 20 = 15 > 10, should not merge
            "*": {1},
        }

        result = merge_index_of_consecutive_indel(mutation_loci)

        # Should merge positions 5-12 for insertions
        expected_insertions = {5, 6, 7, 8, 9, 10, 11, 12}
        assert result["+"] == expected_insertions

        # Should not merge distant deletions
        assert result["-"] == {20, 35}

        # Point mutations unchanged
        assert result["*"] == {1}

    def test_merge_consecutive_already_filled(self):
        """Test when gap is already filled with indels."""
        mutation_loci = {
            "+": {5, 6, 7, 8, 9, 10, 15},  # Consecutive 5-10, then 15
            "-": {1},
            "*": {1},
        }

        result = merge_index_of_consecutive_indel(mutation_loci)

        # Should not add extra positions since 5-10 are already consecutive
        # 15 is too far (15 - 10 = 5 < 10 but gap already filled)
        assert 11 in result["+"]  # Should add positions between 10 and 15
        assert 12 in result["+"]
        assert 13 in result["+"]
        assert 14 in result["+"]

    def test_merge_consecutive_enrichment_both_ends(self):
        """Test enrichment logic when mutations are dense at both ends."""
        # Create a scenario where we have mutations within 10 bases on both sides
        mutation_loci = {
            "+": {1, 2, 3, 4, 5, 10, 15, 16, 17, 18, 19, 20},  # Dense at both ends
            "-": {1},
            "*": {1},
        }

        result = merge_index_of_consecutive_indel(mutation_loci)

        # Should fill the gap between dense regions
        for i in range(6, 10):  # Should fill 6-9
            assert i in result["+"]
        for i in range(11, 15):  # Should fill 11-14
            assert i in result["+"]

    def test_merge_consecutive_no_merging_needed(self):
        """Test when no merging is needed."""
        mutation_loci = {
            "+": {1, 20, 40},  # All far apart
            "-": {5, 25, 45},  # All far apart
            "*": {10},
        }

        result = merge_index_of_consecutive_indel(mutation_loci)

        # Should remain unchanged
        assert result["+"] == {1, 20, 40}
        assert result["-"] == {5, 25, 45}
        assert result["*"] == {10}

    def test_merge_consecutive_empty_sets(self):
        """Test with empty mutation sets."""
        mutation_loci = {"+": set(), "-": set(), "*": set()}

        result = merge_index_of_consecutive_indel(mutation_loci)

        assert result["+"] == set()
        assert result["-"] == set()
        assert result["*"] == set()


class TestSplitKmer:
    """Test split_kmer function."""

    def test_split_kmer_basic(self):
        """Test basic k-mer splitting."""
        indels = {"+": np.array([1, 2, 3, 4, 5]), "-": np.array([5, 4, 3, 2, 1])}
        kmer = 3

        result = split_kmer(indels, kmer)

        assert "+" in result
        assert "-" in result
        assert len(result["+"]) == 5  # Same length as input
        assert len(result["-"]) == 5

    def test_split_kmer_center_positions(self):
        """Test k-mer extraction for center positions."""
        indels = {"+": np.array([1, 2, 3, 4, 5, 6, 7])}
        kmer = 3

        result = split_kmer(indels, kmer)

        # For position 2 (index 2), should get [2, 3, 4]
        expected_center = np.array([2, 3, 4])
        np.testing.assert_array_equal(result["+"][2], expected_center)

    def test_split_kmer_edge_positions(self):
        """Test k-mer extraction for edge positions."""
        indels = {"+": np.array([1, 2, 3, 4, 5])}
        kmer = 3

        result = split_kmer(indels, kmer)

        # For position 0, should get zeros (padding)
        expected_edge = np.array([0, 0, 0])
        np.testing.assert_array_equal(result["+"][0], expected_edge)

    def test_split_kmer_even_size(self):
        """Test k-mer splitting with even k-mer size."""
        indels = {"+": np.array([1, 2, 3, 4, 5, 6])}
        kmer = 4  # Even number

        result = split_kmer(indels, kmer)

        # Should handle even k-mer size correctly
        assert len(result["+"]) == 6
        assert all(len(kmer_array) == 4 for kmer_array in result["+"])

    def test_split_kmer_single_position(self):
        """Test k-mer splitting with single position."""
        indels = {"+": np.array([5])}
        kmer = 3

        result = split_kmer(indels, kmer)

        # Should return zero-filled array for single position
        expected = np.array([0, 0, 0])
        np.testing.assert_array_equal(result["+"][0], expected)


class TestDiscardErrors:
    """Test discard_errors function."""

    def test_discard_errors_basic(self):
        """Test basic error removal."""
        loci = {"+": {1, 2, 3, 4, 5}, "-": {6, 7, 8, 9, 10}, "*": {11, 12, 13, 14, 15}}
        errors = {"+": {2, 4}, "-": {7, 9}, "*": {12, 14}}

        result = discard_errors(loci, errors)

        assert result["+"] == {1, 3, 5}
        assert result["-"] == {6, 8, 10}
        assert result["*"] == {11, 13, 15}

    def test_discard_errors_no_errors(self):
        """Test when no errors need to be removed."""
        loci = {"+": {1, 2, 3}, "-": {4, 5, 6}, "*": {7, 8, 9}}
        errors = {"+": set(), "-": set(), "*": set()}

        result = discard_errors(loci, errors)

        assert result["+"] == {1, 2, 3}
        assert result["-"] == {4, 5, 6}
        assert result["*"] == {7, 8, 9}

    def test_discard_errors_all_errors(self):
        """Test when all loci are errors."""
        loci = {"+": {1, 2}, "-": {3, 4}, "*": {5, 6}}
        errors = {"+": {1, 2}, "-": {3, 4}, "*": {5, 6}}

        result = discard_errors(loci, errors)

        assert result["+"] == set()
        assert result["-"] == set()
        assert result["*"] == set()

    def test_discard_errors_overlapping_errors(self):
        """Test with errors not in loci (should be ignored)."""
        loci = {"+": {1, 3, 5}, "-": {2, 4, 6}, "*": {7, 8, 9}}
        errors = {
            "+": {1, 2, 7},  # 2 and 7 not in loci["+"]
            "-": {3, 4, 8},  # 3 and 8 not in loci["-"]
            "*": {1, 8, 10},  # 1 and 10 not in loci["*"]
        }

        result = discard_errors(loci, errors)

        assert result["+"] == {3, 5}  # Only 1 removed
        assert result["-"] == {2, 6}  # Only 4 removed
        assert result["*"] == {7, 9}  # Only 8 removed


class TestMergeLoci:
    """Test merge_loci function."""

    def test_merge_loci_basic(self):
        """Test basic loci merging."""
        dissimilar_loci = {"+": {1, 2, 3}, "-": {4, 5}, "*": {6}}
        anomal_loci = {"+": {3, 4, 5}, "-": {5, 6, 7}, "*": {7, 8}}

        result = merge_loci(dissimilar_loci, anomal_loci)

        assert result["+"] == {1, 2, 3, 4, 5}
        assert result["-"] == {4, 5, 6, 7}
        assert result["*"] == {6, 7, 8}

    def test_merge_loci_empty_sets(self):
        """Test merging with empty sets."""
        dissimilar_loci = {"+": {1, 2}, "-": set(), "*": {3}}
        anomal_loci = {"+": set(), "-": {4, 5}, "*": set()}

        result = merge_loci(dissimilar_loci, anomal_loci)

        assert result["+"] == {1, 2}
        assert result["-"] == {4, 5}
        assert result["*"] == {3}

    def test_merge_loci_identical_sets(self):
        """Test merging identical sets."""
        loci = {"+": {1, 2, 3}, "-": {4, 5, 6}, "*": {7, 8, 9}}

        result = merge_loci(loci, loci)

        assert result["+"] == {1, 2, 3}
        assert result["-"] == {4, 5, 6}
        assert result["*"] == {7, 8, 9}


class TestAddKnockinLoci:
    """Test add_knockin_loci function."""

    def test_add_knockin_loci_basic(self):
        """Test basic knockin loci addition."""
        candidate_loci = {"+": {1, 2, 3}, "-": {4, 5, 6}, "*": {7, 8, 9}}
        knockin_loci = {10, 11, 12}

        result = add_knockin_loci(candidate_loci, knockin_loci)

        assert result["+"] == {1, 2, 3, 10, 11, 12}
        assert result["-"] == {4, 5, 6, 10, 11, 12}
        assert result["*"] == {7, 8, 9, 10, 11, 12}

    def test_add_knockin_loci_empty_knockin(self):
        """Test with empty knockin loci."""
        candidate_loci = {"+": {1, 2, 3}, "-": {4, 5, 6}, "*": {7, 8, 9}}
        knockin_loci = set()

        result = add_knockin_loci(candidate_loci, knockin_loci)

        assert result["+"] == {1, 2, 3}
        assert result["-"] == {4, 5, 6}
        assert result["*"] == {7, 8, 9}

    def test_add_knockin_loci_overlapping(self):
        """Test with overlapping knockin loci."""
        candidate_loci = {"+": {1, 2, 3}, "-": {2, 3, 4}, "*": {3, 4, 5}}
        knockin_loci = {3, 6, 7}

        result = add_knockin_loci(candidate_loci, knockin_loci)

        assert result["+"] == {1, 2, 3, 6, 7}
        assert result["-"] == {2, 3, 4, 6, 7}
        assert result["*"] == {3, 4, 5, 6, 7}


class TestTransposeMutationLoci:
    """Test transpose_mutation_loci function."""

    def test_transpose_mutation_loci_basic(self):
        """Test basic mutation loci transposition."""
        mutation_loci = {"+": {0, 2, 4}, "-": {1, 3}, "*": {2, 4, 5}}
        sequence = "ACGTAC"  # Length 6

        result = transpose_mutation_loci(mutation_loci, sequence)

        assert len(result) == 6
        assert result[0] == {"+"}
        assert result[1] == {"-"}
        assert result[2] == {"+", "*"}
        assert result[3] == {"-"}
        assert result[4] == {"+", "*"}
        assert result[5] == {"*"}

    def test_transpose_mutation_loci_empty(self):
        """Test transposition with empty mutation loci."""
        mutation_loci = {"+": set(), "-": set(), "*": set()}
        sequence = "ACGT"

        result = transpose_mutation_loci(mutation_loci, sequence)

        assert len(result) == 4
        assert all(loci == set() for loci in result)

    def test_transpose_mutation_loci_single_position(self):
        """Test transposition with single mutation position."""
        mutation_loci = {"+": {2}, "-": set(), "*": set()}
        sequence = "ACGTG"

        result = transpose_mutation_loci(mutation_loci, sequence)

        assert len(result) == 5
        assert result[2] == {"+"}
        assert all(result[i] == set() for i in [0, 1, 3, 4])

    def test_transpose_mutation_loci_out_of_bounds(self):
        """Test transposition with out-of-bounds indices."""
        mutation_loci = {
            "+": {0, 5, 10},  # 5 and 10 are out of bounds for length 4 sequence
            "-": set(),
            "*": set(),
        }
        sequence = "ACGT"  # Length 4

        result = transpose_mutation_loci(mutation_loci, sequence)

        assert len(result) == 4
        assert result[0] == {"+"}
        assert all(result[i] == set() for i in [1, 2, 3])  # Out-of-bounds indices ignored


if __name__ == "__main__":
    pytest.main([__file__])
