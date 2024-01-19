from __future__ import annotations

import gzip
from pathlib import Path

#################################################
# save_concatenated_fastx
#################################################


def extract_extention(file_path: Path) -> str:
    suffixes = file_path.suffixes
    return "".join(suffixes[-2:]) if len(suffixes) >= 2 else suffixes[0]


def is_gzip_file(file_name: Path) -> bool:
    """Check if a file is a GZip compressed file."""
    try:
        with file_name.open("rb") as f:
            return f.read(2) == b"\x1f\x8b"
    except IOError:
        return False


def save_fastq_as_gzip(TEMPDIR: Path, path_fastx: list[Path], barcode: str) -> None:
    """Merge gzip and non-gzip files into a single gzip file."""
    with gzip.open(Path(TEMPDIR, barcode, "fastq", f"{barcode}.fastq.gz"), "wb") as merged_file:
        for file_name in path_fastx:
            if is_gzip_file(file_name):
                with gzip.open(file_name, "rb") as f:
                    merged_file.write(f.read())
            else:
                with open(file_name, "r") as f:
                    merged_file.write(f.read().encode())


def save_concatenated_fastx(TEMPDIR: Path, directory: str) -> None:
    fastx_suffix = {".fa", ".fq", ".fasta", ".fastq", ".fa.gz", ".fq.gz", ".fasta.gz", ".fastq.gz"}
    path_directory = Path(directory)
    barcode = path_directory.stem
    path_fastx = [path for path in path_directory.iterdir() if extract_extention(path) in fastx_suffix]
    save_fastq_as_gzip(TEMPDIR, path_fastx, barcode)
