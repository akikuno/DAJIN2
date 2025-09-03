from __future__ import annotations

from pathlib import Path
from typing import Any

from DAJIN2.core import preprocess
from DAJIN2.utils.bed_handler import bed_to_genome_coordinates


def get_genome_coordinates_from_bed(path_bed: str | Path) -> dict[str, Any]:
    return bed_to_genome_coordinates(path_bed)


def get_genome_coordinates_from_server(
    genome: str, gggenome_url: str, goldenpath_url: str, sequence: str
) -> dict[str, Any]:
    genome_coordinates = {
        "genome": genome,
        "chrom_size": 0,
        "chrom": "control",
        "start": 0,
        "end": len(sequence) - 1,
        "strand": "+",
    }
    # Check if genome coordinates are available (from --genome or --genome-coordinate)
    genome_coordinates = preprocess.fetch_coordinates(genome, gggenome_url, sequence)

    chrom = genome_coordinates["chrom"]
    genome_coordinates["chrom_size"] = preprocess.fetch_chromosome_size(genome, chrom, goldenpath_url)

    return genome_coordinates
