"""
Tests for input formatting and validation in preprocess module.
"""

from collections import defaultdict
from pathlib import Path
from unittest.mock import patch

import pytest

from DAJIN2.core.preprocess.infrastructure.input_formatter import (
    FormattedInputs,
    check_caches,
    convert_input_paths_to_posix,
    create_temporal_directory,
    format_inputs,
    get_genome_coordinates,
    parse_arguments,
)


class TestFormattedInputs:
    """Test FormattedInputs dataclass."""

    def test_formatted_inputs_creation(self):
        """Test basic creation of FormattedInputs."""
        formatted = FormattedInputs(
            path_sample=Path("/path/sample"),
            path_control=Path("/path/control"),
            path_allele=Path("/path/allele.fa"),
            sample_name="test_sample",
            control_name="test_control",
            fasta_alleles={"control": "ACGT"},
            tempdir=Path("/tmp/test"),
            genome_coordinates={"genome": "mm39"},
            threads=4,
            uuid="test-uuid",
        )

        assert formatted.sample_name == "test_sample"
        assert formatted.control_name == "test_control"
        assert formatted.threads == 4
        assert formatted.fasta_alleles == {"control": "ACGT"}

    def test_formatted_inputs_frozen(self):
        """Test that FormattedInputs is frozen (immutable)."""
        formatted = FormattedInputs(
            path_sample=Path("/path/sample"),
            path_control=Path("/path/control"),
            path_allele=Path("/path/allele.fa"),
            sample_name="test_sample",
            control_name="test_control",
            fasta_alleles={"control": "ACGT"},
            tempdir=Path("/tmp/test"),
            genome_coordinates={},
            threads=1,
            uuid="test-uuid",
        )

        # Should not be able to modify frozen dataclass
        with pytest.raises(AttributeError):
            formatted.sample_name = "new_name"


class TestParseArguments:
    """Test parse_arguments function."""

    @patch("uuid.uuid4")
    def test_parse_arguments_basic(self, mock_uuid):
        """Test basic argument parsing."""
        mock_uuid.return_value.hex = "test-uuid-123"

        arguments = {
            "sample": "/path/to/sample",
            "control": "/path/to/control",
            "allele": "/path/to/allele.fa",
            "name": "test_experiment",
            "threads": 4,
        }

        result = parse_arguments(arguments)

        assert len(result) == 7
        path_sample, path_control, path_allele, name, threads, genome_urls, uuid_hex = result

        assert path_sample == Path("/path/to/sample")
        assert path_control == Path("/path/to/control")
        assert path_allele == Path("/path/to/allele.fa")
        assert name == "test_experiment"
        assert threads == 4
        assert uuid_hex == "test-uuid-123"

    @patch("uuid.uuid4")
    def test_parse_arguments_with_genome(self, mock_uuid):
        """Test argument parsing with genome information."""
        mock_uuid.return_value.hex = "test-uuid-456"

        arguments = {
            "sample": "/path/to/sample",
            "control": "/path/to/control",
            "allele": "/path/to/allele.fa",
            "name": "test_experiment",
            "threads": 2,
            "genome": "mm39",
            "blat": "http://blat.url",
            "goldenpath": "http://goldenpath.url",
        }

        result = parse_arguments(arguments)
        path_sample, path_control, path_allele, name, threads, genome_urls, uuid_hex = result

        expected_genome_urls = {"genome": "mm39", "blat": "http://blat.url", "goldenpath": "http://goldenpath.url"}

        assert dict(genome_urls) == expected_genome_urls

    @patch("uuid.uuid4")
    def test_parse_arguments_no_genome(self, mock_uuid):
        """Test argument parsing without genome information."""
        mock_uuid.return_value.hex = "test-uuid-789"

        arguments = {
            "sample": "/path/to/sample",
            "control": "/path/to/control",
            "allele": "/path/to/allele.fa",
            "name": "test_experiment",
            "threads": 1,
        }

        result = parse_arguments(arguments)
        path_sample, path_control, path_allele, name, threads, genome_urls, uuid_hex = result

        # genome_urls should be empty defaultdict
        assert len(genome_urls) == 0
        assert isinstance(genome_urls, defaultdict)


class TestConvertInputPathsToPosix:
    """Test convert_input_paths_to_posix function."""

    @patch("DAJIN2.core.preprocess.input_formatter.io.convert_to_posix")
    def test_convert_input_paths_to_posix(self, mock_convert):
        """Test path conversion to POSIX format."""
        mock_convert.side_effect = lambda x: f"posix_{x}"

        sample, control, allele = convert_input_paths_to_posix(
            "C:\\Windows\\sample", "C:\\Windows\\control", "C:\\Windows\\allele.fa"
        )

        assert sample == "posix_C:\\Windows\\sample"
        assert control == "posix_C:\\Windows\\control"
        assert allele == "posix_C:\\Windows\\allele.fa"

        assert mock_convert.call_count == 3


class TestCreateTemporalDirectory:
    """Test create_temporal_directory function."""

    @patch("DAJIN2.core.preprocess.input_formatter.config.TEMP_ROOT_DIR", "/tmp/dajin")
    @patch("pathlib.Path.mkdir")
    def test_create_temporal_directory(self, mock_mkdir):
        """Test temporal directory creation."""
        result = create_temporal_directory("test_experiment", "control_sample")

        expected_path = Path("/tmp/dajin/test_experiment")
        assert result == expected_path

        # Verify mkdir was called with correct path
        mock_mkdir.assert_called_once_with(parents=True, exist_ok=True)

    @patch("DAJIN2.core.preprocess.input_formatter.config.TEMP_ROOT_DIR", "/custom/temp")
    @patch("pathlib.Path.mkdir")
    def test_create_temporal_directory_custom_root(self, mock_mkdir):
        """Test temporal directory creation with custom root."""
        result = create_temporal_directory("my_exp", "my_control")

        expected_path = Path("/custom/temp/my_exp")
        assert result == expected_path


class TestCheckCaches:
    """Test check_caches function."""

    @patch("DAJIN2.core.preprocess.input_formatter.preprocess.exists_cached_hash")
    @patch("DAJIN2.core.preprocess.input_formatter.preprocess.exists_cached_genome")
    def test_check_caches_both_exist(self, mock_cached_genome, mock_cached_hash):
        """Test when both caches exist."""
        mock_cached_hash.return_value = True
        mock_cached_genome.return_value = True

        tempdir = Path("/tmp/test")
        path_allele = "/path/to/allele.fa"
        genome_url = "mm39"

        result = check_caches(tempdir, path_allele, genome_url)

        assert result is True
        mock_cached_hash.assert_called_once_with(tempdir=tempdir, path=path_allele)
        mock_cached_genome.assert_called_once_with(tempdir=tempdir, genome=genome_url)

    @patch("DAJIN2.core.preprocess.input_formatter.preprocess.exists_cached_hash")
    @patch("DAJIN2.core.preprocess.input_formatter.preprocess.exists_cached_genome")
    def test_check_caches_one_missing(self, mock_cached_genome, mock_cached_hash):
        """Test when one cache is missing."""
        mock_cached_hash.return_value = True
        mock_cached_genome.return_value = False

        tempdir = Path("/tmp/test")
        path_allele = "/path/to/allele.fa"
        genome_url = "hg38"

        result = check_caches(tempdir, path_allele, genome_url)

        assert result is False

    @patch("DAJIN2.core.preprocess.input_formatter.preprocess.exists_cached_hash")
    @patch("DAJIN2.core.preprocess.input_formatter.preprocess.exists_cached_genome")
    def test_check_caches_both_missing(self, mock_cached_genome, mock_cached_hash):
        """Test when both caches are missing."""
        mock_cached_hash.return_value = False
        mock_cached_genome.return_value = False

        tempdir = Path("/tmp/test")
        path_allele = "/path/to/allele.fa"
        genome_url = ""

        result = check_caches(tempdir, path_allele, genome_url)

        assert result is False


class TestGetGenomeCoordinates:
    """Test get_genome_coordinates function."""

    def test_get_genome_coordinates_no_genome(self):
        """Test genome coordinates when no genome is specified."""
        genome_urls = {"genome": ""}
        fasta_alleles = {"control": "ACGTACGTACGT"}
        is_cache_genome = False
        tempdir = Path("/tmp/test")

        result = get_genome_coordinates(genome_urls, fasta_alleles, is_cache_genome, tempdir)

        expected = {
            "genome": "",
            "chrom_size": 0,
            "chrom": "control",
            "start": 0,
            "end": 11,  # len("ACGTACGTACGT") - 1
            "strand": "+",
        }

        assert result == expected

    @patch("DAJIN2.core.preprocess.input_formatter.io.read_jsonl")
    def test_get_genome_coordinates_with_cache(self, mock_read_jsonl):
        """Test genome coordinates with cached data."""
        cached_coords = {
            "genome": "mm39",
            "chrom": "chr1",
            "start": 1000,
            "end": 2000,
            "strand": "+",
            "chrom_size": 195471971,
        }
        mock_read_jsonl.return_value = iter([cached_coords])

        genome_urls = {"genome": "mm39"}
        fasta_alleles = {"control": "ACGT"}
        is_cache_genome = True
        tempdir = Path("/tmp/test")

        result = get_genome_coordinates(genome_urls, fasta_alleles, is_cache_genome, tempdir)

        assert result == cached_coords
        mock_read_jsonl.assert_called_once_with(Path(tempdir, "cache", "genome_coordinates.jsonl"))

    @patch("DAJIN2.core.preprocess.input_formatter.preprocess.fetch_coordinates")
    @patch("DAJIN2.core.preprocess.input_formatter.preprocess.fetch_chromosome_size")
    @patch("DAJIN2.core.preprocess.input_formatter.io.write_jsonl")
    def test_get_genome_coordinates_fetch_new(self, mock_write_jsonl, mock_fetch_size, mock_fetch_coords):
        """Test fetching new genome coordinates."""
        mock_fetch_coords.return_value = {"genome": "hg38", "chrom": "chr2", "start": 500, "end": 1500, "strand": "-"}
        mock_fetch_size.return_value = 242193529

        genome_urls = {"genome": "hg38", "blat": "http://blat.url", "goldenpath": "http://goldenpath.url"}
        fasta_alleles = {"control": "TTTTAAAA"}
        is_cache_genome = False
        tempdir = Path("/tmp/test")

        result = get_genome_coordinates(genome_urls, fasta_alleles, is_cache_genome, tempdir)

        expected = {
            "genome": "hg38",
            "chrom": "chr2",
            "start": 500,
            "end": 1500,
            "strand": "-",
            "chrom_size": 242193529,
        }

        assert result == expected
        mock_fetch_coords.assert_called_once()
        mock_fetch_size.assert_called_once()
        mock_write_jsonl.assert_called_once()


class TestFormatInputs:
    """Test format_inputs integration function."""

    @patch("DAJIN2.core.preprocess.input_formatter.parse_arguments")
    @patch("DAJIN2.core.preprocess.input_formatter.convert_input_paths_to_posix")
    @patch("DAJIN2.core.preprocess.input_formatter.fastx_handler.extract_filename")
    @patch("DAJIN2.core.preprocess.input_formatter.fastx_handler.dictionize_allele")
    @patch("DAJIN2.core.preprocess.input_formatter.create_temporal_directory")
    @patch("DAJIN2.core.preprocess.input_formatter.check_caches")
    @patch("DAJIN2.core.preprocess.input_formatter.get_genome_coordinates")
    @patch("DAJIN2.core.preprocess.input_formatter.io.sanitize_name")
    def test_format_inputs_integration(
        self,
        mock_sanitize,
        mock_get_coords,
        mock_check_caches,
        mock_create_tempdir,
        mock_dictionize,
        mock_extract_filename,
        mock_convert_paths,
        mock_parse_args,
    ):
        """Test complete format_inputs integration."""
        # Setup mocks
        mock_parse_args.return_value = (
            Path("/path/sample"),
            Path("/path/control"),
            Path("/path/allele.fa"),
            "test_exp",
            4,
            {"genome": ""},
            "uuid-123",
        )
        mock_convert_paths.return_value = ("/posix/sample", "/posix/control", "/posix/allele.fa")
        mock_extract_filename.side_effect = ["sample_name", "control_name"]
        mock_dictionize.return_value = {"control": "ACGT", "target": "TTTT"}
        mock_create_tempdir.return_value = Path("/tmp/test_exp")
        mock_check_caches.return_value = False
        mock_get_coords.return_value = {"genome": "", "chrom": "control"}
        mock_sanitize.side_effect = lambda x: f"clean_{x}"

        arguments = {
            "sample": "/path/sample",
            "control": "/path/control",
            "allele": "/path/allele.fa",
            "name": "test_exp",
            "threads": 4,
        }

        result = format_inputs(arguments)

        # Verify result
        assert isinstance(result, FormattedInputs)
        assert result.sample_name == "clean_sample_name"
        assert result.control_name == "clean_control_name"
        assert result.threads == 4
        assert result.uuid == "uuid-123"
        assert result.fasta_alleles == {"control": "ACGT", "target": "TTTT"}

        # Verify all functions were called
        mock_parse_args.assert_called_once_with(arguments)
        mock_convert_paths.assert_called_once()
        assert mock_extract_filename.call_count == 2
        mock_dictionize.assert_called_once()
        mock_create_tempdir.assert_called_once()
        mock_check_caches.assert_called_once()
        mock_get_coords.assert_called_once()
        assert mock_sanitize.call_count == 2


if __name__ == "__main__":
    pytest.main([__file__])
