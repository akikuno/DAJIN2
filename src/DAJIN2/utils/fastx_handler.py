from __future__ import annotations

import re
import gzip
from pathlib import Path

import mappy


#################################################
# Helper function
#################################################


def sanitize_filename(path_file: Path | str) -> str:
    """
    Sanitize the path_file by replacing invalid characters on Windows OS with '-'
    """
    path_file = str(path_file).lstrip()
    if not path_file:
        raise ValueError("Provided FASTA/FASTQ is empty or consists only of whitespace")
    return re.sub(r'[\\/:?.,\'"<>| ]', "-", path_file)


#################################################
# Extract filename
#################################################


def extract_filename(path_fasta: Path | str) -> str:
    filename = Path(path_fasta).name
    filename = re.sub(r"\..*$", "", filename)  # Remove file extension
    return sanitize_filename(filename)


#################################################
# Convert allele file to dictionary type fasta format
#################################################


def dictionize_allele(path_fasta: str | Path) -> dict[str, str]:
    return {sanitize_filename(name): seq.upper() for name, seq, _ in mappy.fastx_read(str(path_fasta))}


#################################################
# Export fasta files as single-FASTA format
#################################################


def export_fasta_files(TEMPDIR: Path, FASTA_ALLELES: dict, NAME: str) -> None:
    """+ Save multiple FASTAs in separate single-FASTA format files."""
    for identifier, sequence in FASTA_ALLELES.items():
        contents = "\n".join([">" + identifier, sequence]) + "\n"
        output_fasta = Path(TEMPDIR, NAME, "fasta", f"{identifier}.fasta")
        output_fasta.write_text(contents)


#################################################
# save_concatenated_fastx
#################################################


def extract_extention(path_file: Path) -> str:
    suffixes = path_file.suffixes
    return "".join(suffixes)


def is_gzip_file(path_file: Path) -> bool:
    """Check if a file is a GZip compressed file."""
    try:
        with path_file.open("rb") as f:
            return f.read(2) == b"\x1f\x8b"
    except IOError:
        return False


def save_fastq_as_gzip(TEMPDIR: Path, path_fastx: list[Path], barcode: str) -> None:
    """Merge gzip and non-gzip files into a single gzip file."""
    with gzip.open(Path(TEMPDIR, barcode, "fastq", f"{barcode}.fastq.gz"), "wb") as merged_file:
        for path_file in path_fastx:
            if is_gzip_file(path_file):
                with gzip.open(path_file, "rb") as f:
                    merged_file.write(f.read())
            else:
                with open(path_file, "r") as f:
                    merged_file.write(f.read().encode())


def save_concatenated_fastx(TEMPDIR: Path, directory: str) -> None:
    fastx_suffix = {".fa", ".fq", ".fasta", ".fastq", ".fa.gz", ".fq.gz", ".fasta.gz", ".fastq.gz"}
    path_directory = Path(directory)
    barcode = path_directory.stem
    path_fastx = [path for path in path_directory.iterdir() if extract_extention(path) in fastx_suffix]
    save_fastq_as_gzip(TEMPDIR, path_fastx, barcode)
