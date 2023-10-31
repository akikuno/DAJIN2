from __future__ import annotations

import re
from pathlib import Path

import mappy

########################################################################
# Helper function
########################################################################


def _sanitize_name(name: str) -> str:
    """
    Sanitize the name by replacing invalid characters with '-'
    """
    name = name.lstrip()
    if not name:
        raise ValueError("Provided FASTA/FASTQ is empty or consists only of whitespace")
    return re.sub(r'[\\/:?.,\'"<>| ]', "-", name)


########################################################################
# Extract basename
########################################################################


def extract_basename(fastq_path: str) -> str:
    name = Path(fastq_path).name
    name = re.sub(r"\..*$", "", name)  # Remove file extension
    return _sanitize_name(name)


########################################################################
# Convert allele file to dictionary type fasta format
########################################################################


def dictionize_allele(path_fasta: str | Path) -> dict[str, str]:
    return {_sanitize_name(name): seq.upper() for name, seq, _ in mappy.fastx_read(str(path_fasta))}


########################################################################
# Export fasta files as single-FASTA format
########################################################################


def export_fasta_files(TEMPDIR: Path, FASTA_ALLELES: dict, NAME: str) -> None:
    """
    This function exports FASTA files in single-FASTA format.

    :param TEMPDIR: Temporary directory Path object where the output files will be saved.
    :param FASTA_ALLELES: Dictionary containing identifier and sequence pairs.
    :param NAME: Name to be included in the output path.
    """
    for identifier, sequence in FASTA_ALLELES.items():
        contents = "\n".join([">" + identifier, sequence]) + "\n"
        output_fasta = Path(TEMPDIR, NAME, "fasta", f"{identifier}.fasta")
        output_fasta.write_text(contents)
