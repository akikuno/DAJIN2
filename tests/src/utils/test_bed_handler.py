"""
Tests for BED file handling utilities.
"""

import tempfile
from pathlib import Path
from unittest.mock import patch

import pytest

from DAJIN2.utils.bed_handler import (
    BEDError,
    bed_to_genome_coordinates,
    format_coordinates_summary,
    parse_bed_file,
    validate_bed_coordinates,
)


class TestParseBedFile:
    """Test parse_bed_file function."""

    def test_parse_bed_file_basic(self):
        """Test basic BED file parsing."""
        bed_content = "chr1\t100\t200\t248956422\t0\t+"

        with tempfile.NamedTemporaryFile(mode="w", suffix=".bed", delete=False) as f:
            f.write(bed_content)
            f.flush()

            result = parse_bed_file(f.name)

        Path(f.name).unlink()  # Clean up

        expected = [
            {"chrom": "chr1", "start": 100, "end": 200, "name": "248956422", "chrom_size": 248956422, "strand": "+"}
        ]

        assert result == expected

    def test_parse_bed_file_minimal(self):
        """Test BED file with minimal 6 columns (BED6 format)."""
        bed_content = "chr2\t1000\t2000\t195471971\t0\t+"

        with tempfile.NamedTemporaryFile(mode="w", suffix=".bed", delete=False) as f:
            f.write(bed_content)
            f.flush()

            result = parse_bed_file(f.name)

        Path(f.name).unlink()

        expected = [
            {"chrom": "chr2", "start": 1000, "end": 2000, "name": "195471971", "chrom_size": 195471971, "strand": "+"}
        ]

        assert result == expected

    def test_parse_bed_file_multiple_intervals(self):
        """Test BED file with multiple intervals."""
        bed_content = (
            "chr1\t100\t200\t248956422\t0\t+\nchr2\t300\t400\t195471971\t0\t-\nchr3\t500\t600\t159345973\t0\t+"
        )

        with tempfile.NamedTemporaryFile(mode="w", suffix=".bed", delete=False) as f:
            f.write(bed_content)
            f.flush()

            result = parse_bed_file(f.name)

        Path(f.name).unlink()

        assert len(result) == 3
        assert result[0]["chrom"] == "chr1"
        assert result[1]["chrom"] == "chr2"
        assert result[2]["chrom"] == "chr3"

    def test_parse_bed_file_with_comments(self):
        """Test BED file with comments and track lines."""
        bed_content = """# Comment line
track name="test track"
chr1\t100\t200\t248956422\t0\t+
#another comment
chr2\t300\t400\t195471971\t0\t-"""

        with tempfile.NamedTemporaryFile(mode="w", suffix=".bed", delete=False) as f:
            f.write(bed_content)
            f.flush()

            result = parse_bed_file(f.name)

        Path(f.name).unlink()

        assert len(result) == 2
        assert result[0]["chrom"] == "chr1"
        assert result[1]["chrom"] == "chr2"

    def test_parse_bed_file_invalid_coordinates(self):
        """Test BED file with invalid coordinates."""
        bed_content = "chr1\t200\t100\t248956422\t0\t+"  # end < start

        with tempfile.NamedTemporaryFile(mode="w", suffix=".bed", delete=False) as f:
            f.write(bed_content)
            f.flush()

            with pytest.raises(BEDError, match="Invalid end position"):
                parse_bed_file(f.name)

        Path(f.name).unlink()

    def test_parse_bed_file_invalid_format_too_few_columns(self):
        """Test BED file with too few columns."""
        bed_content = "chr1\t100"  # Only 2 columns

        with tempfile.NamedTemporaryFile(mode="w", suffix=".bed", delete=False) as f:
            f.write(bed_content)
            f.flush()

            with pytest.raises(BEDError, match="DAJIN2 requires BED6 format"):
                parse_bed_file(f.name)

        Path(f.name).unlink()

    def test_parse_bed_file_bed3_format_rejected(self):
        """Test that BED3 format is rejected (requires strand)."""
        bed_content = "chr1\t100\t200"  # BED3 format without strand

        with tempfile.NamedTemporaryFile(mode="w", suffix=".bed", delete=False) as f:
            f.write(bed_content)
            f.flush()

            with pytest.raises(BEDError, match="DAJIN2 requires BED6 format"):
                parse_bed_file(f.name)

        Path(f.name).unlink()

    def test_parse_bed_file_bed5_format_rejected(self):
        """Test that BED5 format is rejected (missing strand)."""
        bed_content = "chr1\t100\t200\tfeature\t0"  # BED5 format without strand

        with tempfile.NamedTemporaryFile(mode="w", suffix=".bed", delete=False) as f:
            f.write(bed_content)
            f.flush()

            with pytest.raises(BEDError, match="DAJIN2 requires BED6 format"):
                parse_bed_file(f.name)

        Path(f.name).unlink()

    def test_parse_bed_file_invalid_strand(self):
        """Test BED file with invalid strand."""
        bed_content = "chr1\t100\t200\t248956422\t0\tx"  # Invalid strand 'x'

        with tempfile.NamedTemporaryFile(mode="w", suffix=".bed", delete=False) as f:
            f.write(bed_content)
            f.flush()

            with pytest.raises(BEDError, match="Invalid or missing strand"):
                parse_bed_file(f.name)

        Path(f.name).unlink()

    def test_parse_bed_file_empty_strand(self):
        """Test BED file with empty strand field."""
        bed_content = "chr1\t100\t200\t248956422\t0\t"  # Empty strand

        with tempfile.NamedTemporaryFile(mode="w", suffix=".bed", delete=False) as f:
            f.write(bed_content)
            f.flush()

            with pytest.raises(BEDError, match="Invalid BED format.*Expected 6 fields.*got 5"):
                parse_bed_file(f.name)

        Path(f.name).unlink()

    def test_parse_bed_file_invalid_chromosome_size_string(self):
        """Test BED file with non-integer chromosome size."""
        bed_content = "chr1\t100\t200\tfeature_name\t0\t+"  # String instead of integer

        with tempfile.NamedTemporaryFile(mode="w", suffix=".bed", delete=False) as f:
            f.write(bed_content)
            f.flush()

            with pytest.raises(BEDError, match="Invalid chromosome size format"):
                parse_bed_file(f.name)

        Path(f.name).unlink()

    def test_parse_bed_file_negative_chromosome_size(self):
        """Test BED file with negative chromosome size."""
        bed_content = "chr1\t100\t200\t-123456\t0\t+"  # Negative chromosome size

        with tempfile.NamedTemporaryFile(mode="w", suffix=".bed", delete=False) as f:
            f.write(bed_content)
            f.flush()

            with pytest.raises(BEDError, match="Invalid chromosome size.*must be a positive integer"):
                parse_bed_file(f.name)

        Path(f.name).unlink()

    def test_parse_bed_file_zero_chromosome_size(self):
        """Test BED file with zero chromosome size."""
        bed_content = "chr1\t100\t200\t0\t0\t+"  # Zero chromosome size

        with tempfile.NamedTemporaryFile(mode="w", suffix=".bed", delete=False) as f:
            f.write(bed_content)
            f.flush()

            with pytest.raises(BEDError, match="Invalid chromosome size.*must be a positive integer"):
                parse_bed_file(f.name)

        Path(f.name).unlink()

    def test_parse_bed_file_valid_chromosome_size(self):
        """Test BED file with valid chromosome size."""
        bed_content = "chr1\t100\t200\t248956422\t0\t+"  # Valid chromosome size

        with tempfile.NamedTemporaryFile(mode="w", suffix=".bed", delete=False) as f:
            f.write(bed_content)
            f.flush()

            result = parse_bed_file(f.name)

        Path(f.name).unlink()

        expected = [
            {"chrom": "chr1", "start": 100, "end": 200, "name": "248956422", "chrom_size": 248956422, "strand": "+"}
        ]

        assert result == expected

    def test_parse_bed_file_nonexistent(self):
        """Test parsing non-existent BED file."""
        with pytest.raises(FileNotFoundError):
            parse_bed_file("/nonexistent/file.bed")

    def test_parse_bed_file_empty(self):
        """Test parsing empty BED file."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".bed", delete=False) as f:
            f.write("")  # Empty file
            f.flush()

            with pytest.raises(BEDError, match="No valid intervals found"):
                parse_bed_file(f.name)

        Path(f.name).unlink()


class TestValidateBedCoordinates:
    """Test validate_bed_coordinates function."""

    def test_validate_bed_coordinates_valid(self):
        """Test validation of valid coordinates."""
        intervals = [{"chrom": "chr1", "start": 100, "end": 200}, {"chrom": "chr2", "start": 300, "end": 400}]

        # Should not raise any exception
        validate_bed_coordinates(intervals)

    def test_validate_bed_coordinates_invalid_types(self):
        """Test validation with invalid coordinate types."""
        intervals = [
            {"chrom": "chr1", "start": "100", "end": 200}  # start is string
        ]

        with pytest.raises(BEDError, match="coordinates must be integers"):
            validate_bed_coordinates(intervals)

    def test_validate_bed_coordinates_invalid_order(self):
        """Test validation with start >= end."""
        intervals = [
            {"chrom": "chr1", "start": 200, "end": 200}  # start == end
        ]

        with pytest.raises(BEDError, match="start .* must be less than end"):
            validate_bed_coordinates(intervals)

    def test_validate_bed_coordinates_empty(self):
        """Test validation with empty intervals."""
        with pytest.raises(BEDError, match="No intervals provided"):
            validate_bed_coordinates([])


class TestBedToGenomeCoordinates:
    """Test bed_to_genome_coordinates function."""

    def test_bed_to_genome_coordinates_basic(self):
        """Test conversion of BED to genome coordinates."""
        bed_content = "chr1\t100\t200\t248956422\t0\t+"

        with tempfile.NamedTemporaryFile(mode="w", suffix=".bed", delete=False) as f:
            f.write(bed_content)
            f.flush()

            result = bed_to_genome_coordinates(f.name, "hg38")

        Path(f.name).unlink()

        expected = {
            "genome": "hg38",
            "chrom": "chr1",
            "start": 100,
            "end": 200,
            "strand": "+",
            "chrom_size": 248956422,
        }

        assert result == expected

    def test_bed_to_genome_coordinates_no_genome(self):
        """Test conversion without genome specification."""
        bed_content = "chr2\t500\t600\t195471971\t0\t-"

        with tempfile.NamedTemporaryFile(mode="w", suffix=".bed", delete=False) as f:
            f.write(bed_content)
            f.flush()

            result = bed_to_genome_coordinates(f.name)

        Path(f.name).unlink()

        assert result["genome"] == ""
        assert result["chrom"] == "chr2"
        assert result["start"] == 500
        assert result["end"] == 600
        assert result["strand"] == "-"
        assert result["chrom_size"] == 195471971

    @patch("logging.getLogger")
    def test_bed_to_genome_coordinates_multiple_intervals_warning(self, mock_logger):
        """Test that multiple intervals generate a warning."""
        bed_content = "chr1\t100\t200\t248956422\t0\t+\nchr2\t300\t400\t195471971\t0\t-"

        with tempfile.NamedTemporaryFile(mode="w", suffix=".bed", delete=False) as f:
            f.write(bed_content)
            f.flush()

            result = bed_to_genome_coordinates(f.name, "mm39")

        Path(f.name).unlink()

        # Should use first interval
        assert result["chrom"] == "chr1"
        assert result["start"] == 100
        assert result["end"] == 200

        # Should have logged a warning
        mock_logger.return_value.warning.assert_called_once()


class TestFormatCoordinatesSummary:
    """Test format_coordinates_summary function."""

    def test_format_coordinates_summary_complete(self):
        """Test formatting with complete coordinates."""
        genome_coordinates = {"genome": "hg38", "chrom": "chr1", "start": 100, "end": 200, "strand": "+"}

        result = format_coordinates_summary(genome_coordinates)

        assert result == "hg38 chr1:100-200(+)"

    def test_format_coordinates_summary_no_genome(self):
        """Test formatting without genome."""
        genome_coordinates = {"chrom": "chr2", "start": 500, "end": 600, "strand": "-"}

        result = format_coordinates_summary(genome_coordinates)

        assert result == "chr2:500-600(-)"

    def test_format_coordinates_summary_defaults(self):
        """Test formatting with missing fields."""
        genome_coordinates = {}

        result = format_coordinates_summary(genome_coordinates)

        assert result == "unknown:0-0(+)"
