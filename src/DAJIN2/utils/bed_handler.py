"""
BED file handling utilities for DAJIN2.

This module provides functions to parse BED files and extract genomic coordinates
for use in DAJIN2 analysis.
"""

from __future__ import annotations

from pathlib import Path


class BEDError(Exception):
    """Exception raised for BED file parsing errors."""

    def __init__(self, message: str) -> None:
        super().__init__(message)
        self.message = message


def parse_bed_file(bed_path: str | Path) -> list[dict[str, any]]:
    """
    Parse a BED file and return a list of genomic intervals.

    Args:
        bed_path: Path to the BED file

    Returns:
        List of dictionaries containing genomic intervals with keys:
        - chrom: chromosome name
        - start: 0-based start position (int)
        - end: 0-based end position (int)
        - name: feature name (optional)
        - strand: strand information (optional, + or -)

    Raises:
        BEDError: If BED file is invalid or cannot be parsed
        FileNotFoundError: If BED file does not exist
    """
    bed_path = Path(bed_path)

    if not bed_path.exists():
        raise FileNotFoundError(f"BED file not found: {bed_path}")

    intervals = []

    try:
        with open(bed_path) as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()

                # Skip empty lines and comments
                if not line or line.startswith("#") or line.startswith("track"):
                    continue

                fields = line.split("\t")

                # BED format requires at least 6 fields for DAJIN2: chrom, start, end, name, score, strand
                if len(fields) < 6:
                    raise BEDError(
                        f"Invalid BED format at line {line_num}: "
                        f"DAJIN2 requires BED6 format with strand information. "
                        f"Expected 6 fields (chrom, start, end, name, score, strand), got {len(fields)}. "
                        f"Please use BED6 format with strand column (+/-) for proper sequence orientation."
                    )

                try:
                    chrom = fields[0]
                    start = int(fields[1])  # BED uses 0-based start
                    end = int(fields[2])  # BED uses 1-based end (exclusive)

                    # Validate coordinates
                    if start < 0:
                        raise BEDError(f"Invalid start position at line {line_num}: {start} (must be >= 0)")
                    if end <= start:
                        raise BEDError(f"Invalid end position at line {line_num}: {end} (must be > start)")

                except ValueError as e:
                    raise BEDError(f"Invalid coordinate format at line {line_num}: {e}")

                # Create interval dictionary
                interval = {"chrom": chrom, "start": start, "end": end}

                # Process required BED6 fields - chromosome size validation
                if len(fields) >= 4 and fields[3]:
                    try:
                        chrom_size = int(fields[3])
                        if chrom_size <= 0:
                            raise BEDError(
                                f"Invalid chromosome size at line {line_num}: {chrom_size} "
                                f"(must be a positive integer). Please use the actual chromosome size."
                            )
                        interval["name"] = fields[3]
                        interval["chrom_size"] = chrom_size
                    except ValueError:
                        raise BEDError(
                            f"Invalid chromosome size format at line {line_num}: '{fields[3]}' "
                            f"(must be an integer). Please enter the chromosome size in the 4th column of BED6 format. "
                            f"Example: chr1\t1000000\t1001000\t248956422\t0\t+"
                        )

                # Strand is required (field 6)
                if fields[5] in ["+", "-"]:
                    interval["strand"] = fields[5]
                else:
                    raise BEDError(
                        f"Invalid or missing strand at line {line_num}: '{fields[5]}' "
                        f"(must be '+' or '-'). DAJIN2 requires strand information for proper sequence orientation."
                    )

                intervals.append(interval)

    except Exception as e:
        if isinstance(e, BEDError):
            raise
        else:
            raise BEDError(f"Error reading BED file {bed_path}: {e}")

    if not intervals:
        raise BEDError(f"No valid intervals found in BED file: {bed_path}")

    # Validate coordinates (integrated from validate_bed_coordinates)
    for i, interval in enumerate(intervals):
        start = interval["start"]
        end = interval["end"]

        # Basic validation (already done above, but keeping for safety)
        if not isinstance(start, int) or not isinstance(end, int):
            raise BEDError(f"Interval {i + 1}: coordinates must be integers")

        if start >= end:
            raise BEDError(f"Interval {i + 1}: start ({start}) must be less than end ({end})")

    return intervals


def validate_bed_coordinates(intervals: list[dict[str, any]]) -> None:
    """
    Validate BED coordinates for correctness.

    Args:
        intervals: List of interval dictionaries with 'start' and 'end' keys

    Raises:
        BEDError: If coordinates are invalid or intervals list is empty
    """
    if not intervals:
        raise BEDError("No intervals provided for validation")

    for i, interval in enumerate(intervals):
        start = interval.get("start")
        end = interval.get("end")

        # Check if coordinates are integers
        if not isinstance(start, int) or not isinstance(end, int):
            raise BEDError(f"Interval {i + 1}: coordinates must be integers")

        # Check if start < end
        if start >= end:
            raise BEDError(f"Interval {i + 1}: start ({start}) must be less than end ({end})")


def bed_to_genome_coordinates(bed_path: str | Path, genome: str = "") -> dict[str, any]:
    """
    Convert BED file to DAJIN2 genome_coordinates format.

    For multi-interval BED files, uses the first interval.

    Args:
        bed_path: Path to BED file
        genome: Genome assembly name (e.g., "hg38", "mm39")

    Returns:
        Dictionary in DAJIN2 genome_coordinates format:
        {
            "genome": genome,
            "chrom": chromosome,
            "start": start_position,
            "end": end_position,
            "strand": strand,
            "chrom_size": 0
        }

    Raises:
        BEDError: If BED file is invalid
    """
    intervals = parse_bed_file(bed_path)

    # Use first interval for single-region analysis
    first_interval = intervals[0]

    # Warn if multiple intervals
    if len(intervals) > 1:
        import logging

        logger = logging.getLogger(__name__)
        logger.warning(
            f"BED file contains {len(intervals)} intervals. Using first interval: "
            f"{first_interval['chrom']}:{first_interval['start']}-{first_interval['end']}"
        )

    # Convert to DAJIN2 format
    genome_coordinates = {
        "genome": genome,
        "chrom": first_interval["chrom"],
        "start": first_interval["start"],  # Already 0-based
        "end": first_interval["end"],  # Already exclusive
        "strand": first_interval.get("strand", "+"),
        "chrom_size": first_interval.get("chrom_size", 0),  # Use chromosome size from BED 4th column
    }

    return genome_coordinates


def format_coordinates_summary(genome_coordinates: dict[str, any]) -> str:
    """
    Format genome coordinates for logging/display.

    Args:
        genome_coordinates: Dictionary in DAJIN2 format

    Returns:
        Formatted string representation
    """
    chrom = genome_coordinates.get("chrom", "unknown")
    start = genome_coordinates.get("start", 0)
    end = genome_coordinates.get("end", 0)
    strand = genome_coordinates.get("strand", "+")
    genome = genome_coordinates.get("genome", "")

    coord_str = f"{chrom}:{start}-{end}({strand})"
    if genome:
        coord_str = f"{genome} {coord_str}"

    return coord_str
