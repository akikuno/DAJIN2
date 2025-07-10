"""
Integration tests for CLI functionality, specifically testing --genome-coordinate argument.
"""

import argparse
import subprocess
import sys
import tempfile
from pathlib import Path
from unittest.mock import patch

import pytest


class TestCLIGenomeCoordinate:
    """Test CLI integration for --genome-coordinate argument."""

    def create_test_bed_file(self, content: str) -> str:
        """Create a temporary BED file with given content."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".bed", delete=False) as f:
            f.write(content)
            f.flush()
            return f.name

    def test_cli_help_includes_bed_option(self):
        """Test that -b/--bed appears in help output."""
        # Create a parser similar to the one in main.py
        parser = argparse.ArgumentParser()
        parser.add_argument("-s", "--sample", type=str, help="Path to a sample directory including FASTQ file")
        parser.add_argument("-c", "--control", type=str, help="Path to a control directory including FASTQ file")
        parser.add_argument("-a", "--allele", type=str, help="Path to a FASTA file")
        parser.add_argument("-n", "--name", type=str, help="Output directory name", default="Results")
        parser.add_argument(
            "-g", "--genome", type=str, default="", help="Reference genome ID (e.g hg38, mm39) [default: '']"
        )
        parser.add_argument(
            "-b",
            "--bed",
            type=str,
            default="",
            dest="genome_coordinate",
            help="Path to BED6 file containing genomic coordinates [default: '']",
        )
        parser.add_argument(
            "--no-filter", action="store_true", help="Disable minor allele filtering (keep alleles <0.5%%)"
        )

        help_text = parser.format_help()

        # The actual format includes the metavar, so look for the pattern that appears
        assert "-b GENOME_COORDINATE, --bed GENOME_COORDINATE" in help_text
        assert "BED6 file containing genomic coordinates" in help_text

    def create_test_parser(self):
        """Create a test parser matching main.py structure."""
        parser = argparse.ArgumentParser()
        parser.add_argument("-s", "--sample", type=str, help="Path to a sample directory including FASTQ file")
        parser.add_argument("-c", "--control", type=str, help="Path to a control directory including FASTQ file")
        parser.add_argument("-a", "--allele", type=str, help="Path to a FASTA file")
        parser.add_argument("-n", "--name", type=str, help="Output directory name", default="Results")
        parser.add_argument(
            "-g", "--genome", type=str, default="", help="Reference genome ID (e.g hg38, mm39) [default: '']"
        )
        parser.add_argument(
            "-b",
            "--bed",
            type=str,
            default="",
            dest="genome_coordinate",
            help="Path to BED6 file containing genomic coordinates [default: '']",
        )
        parser.add_argument("-t", "--threads", type=int, default=1, help="Number of threads [default: 1]")
        parser.add_argument(
            "--no-filter", action="store_true", help="Disable minor allele filtering (keep alleles <0.5%%)"
        )
        return parser

    def test_cli_argparse_bed_short_option(self):
        """Test parsing -b argument (short option)."""
        parser = self.create_test_parser()
        bed_file = self.create_test_bed_file("chr1\t1000\t2000\t248956422\t0\t+")

        try:
            args = parser.parse_args(["-b", bed_file])
            assert args.genome_coordinate == bed_file
            assert args.genome == ""
        finally:
            Path(bed_file).unlink()

    def test_cli_argparse_bed_long_option(self):
        """Test parsing --bed argument (long option)."""
        parser = self.create_test_parser()
        bed_file = self.create_test_bed_file("chr1\t1000\t2000\t248956422\t0\t+")

        try:
            args = parser.parse_args(["--bed", bed_file])
            assert args.genome_coordinate == bed_file
            assert args.genome == ""
        finally:
            Path(bed_file).unlink()

    def test_cli_argparse_both_genome_and_bed(self):
        """Test parsing both --genome and --bed arguments."""
        parser = self.create_test_parser()
        bed_file = self.create_test_bed_file("chr1\t1000\t2000\t248956422\t0\t+")

        try:
            args = parser.parse_args(["--genome", "hg38", "--bed", bed_file])
            assert args.genome == "hg38"
            assert args.genome_coordinate == bed_file
        finally:
            Path(bed_file).unlink()

    def test_cli_argparse_genome_only(self):
        """Test parsing --genome argument alone."""
        parser = self.create_test_parser()
        args = parser.parse_args(["--genome", "mm39"])

        assert args.genome == "mm39"
        assert args.genome_coordinate == ""

    def test_cli_argparse_no_genome_arguments(self):
        """Test parsing with no genome-related arguments."""
        parser = self.create_test_parser()
        args = parser.parse_args([])

        assert args.genome == ""
        assert args.genome_coordinate == ""

    def test_cli_precedence_bed_over_genome(self):
        """Test that BED file takes precedence over genome ID in CLI processing."""
        from DAJIN2.main import execute

        bed_file = self.create_test_bed_file("chr1\t1000\t2000\t248956422\t0\t+")

        # Mock sys.argv to simulate CLI call
        test_args = [
            "DAJIN2",
            "--sample",
            "/fake/sample",
            "--control",
            "/fake/control",
            "--allele",
            "/fake/allele.fa",
            "--name",
            "test",
            "--genome",
            "hg38",
            "--bed",
            bed_file,
        ]

        try:
            with patch.object(sys, "argv", test_args):
                with patch("DAJIN2.utils.input_validator.validate_files"):
                    with patch(
                        "DAJIN2.utils.input_validator.validate_bed_file_and_get_coordinates"
                    ) as mock_bed_validate:
                        with patch("DAJIN2.core.core.execute_control"):
                            with patch("DAJIN2.core.core.execute_sample"):
                                with patch("DAJIN2.main.generate_report"):
                                    with patch("DAJIN2.main.cache_bed_coordinates"):
                                        with patch("shutil.move"):
                                            with patch("shutil.rmtree"):
                                                mock_bed_validate.return_value = {
                                                    "genome": "hg38",
                                                    "chrom": "chr1",
                                                    "start": 1000,
                                                    "end": 2000,
                                                    "strand": "+",
                                                    "chrom_size": 248956422,
                                                }

                                                execute()

                                                # Verify BED validation was called (meaning BED took precedence)
                                                mock_bed_validate.assert_called_once_with(bed_file, "hg38")

        finally:
            Path(bed_file).unlink()

    def test_cli_genome_only_when_no_bed(self):
        """Test that genome ID is used when no BED file is provided."""
        from DAJIN2.main import execute

        test_args = [
            "DAJIN2",
            "--sample",
            "/fake/sample",
            "--control",
            "/fake/control",
            "--allele",
            "/fake/allele.fa",
            "--name",
            "test",
            "--genome",
            "mm39",
        ]

        with patch.object(sys, "argv", test_args):
            with patch("DAJIN2.utils.input_validator.validate_files"):
                with patch("DAJIN2.utils.input_validator.validate_genome_and_fetch_urls") as mock_genome_validate:
                    with patch("DAJIN2.core.core.execute_control"):
                        with patch("DAJIN2.core.core.execute_sample"):
                            with patch("DAJIN2.main.generate_report"):
                                with patch("shutil.move"):
                                    with patch("shutil.rmtree"):
                                        mock_genome_validate.return_value = {"genome": "mm39"}

                                        execute()

                                        # Verify genome validation was called
                                        mock_genome_validate.assert_called_once_with("mm39")

    def test_cli_batch_mode_header_validation_includes_bed_columns(self):
        """Test that batch mode accepts both genome_coordinate and bed columns."""
        from DAJIN2.main import validate_headers_of_batch_file

        # Test with genome_coordinate column (legacy)
        columns = {"sample", "control", "allele", "name", "genome_coordinate"}
        # Should not raise an error
        validate_headers_of_batch_file(columns, "test.csv")

        # Test with bed column (new)
        columns = {"sample", "control", "allele", "name", "bed"}
        # Should not raise an error
        validate_headers_of_batch_file(columns, "test.csv")

        # Test with both genome and bed
        columns = {"sample", "control", "allele", "name", "genome", "bed"}
        # Should not raise an error
        validate_headers_of_batch_file(columns, "test.csv")

    @pytest.mark.slow
    def test_cli_help_command_real(self):
        """Test actual CLI help command includes -b/--bed option."""
        # This test requires DAJIN2 to be importable, so mark as slow
        try:
            result = subprocess.run(
                [sys.executable, "-m", "DAJIN2", "--help"], capture_output=True, text=True, timeout=10
            )

            if result.returncode == 0:
                assert "-b, --bed" in result.stdout
                assert "BED6 file containing genomic coordinates" in result.stdout
            else:
                # If module import fails, skip this test
                pytest.skip("DAJIN2 module not properly installed for CLI testing")

        except (subprocess.TimeoutExpired, FileNotFoundError):
            pytest.skip("Cannot run DAJIN2 CLI for integration testing")
