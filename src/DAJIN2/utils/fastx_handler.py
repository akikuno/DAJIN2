from __future__ import annotations

import gzip
import random
import re
import uuid
from pathlib import Path

import mappy
import pysam
from DAJIN2.utils.io import sanitize_name

#################################################
# Extract filename
#################################################


def extract_filename(path_fasta: Path | str) -> str:
    filename = Path(path_fasta).name
    filename = re.sub(r"\..*$", "", filename)  # Remove file extension
    return sanitize_name(filename)


#################################################
# Convert allele file to dictionary type fasta format
#################################################


def dictionize_allele(path_fasta: str | Path) -> dict[str, str]:
    return {sanitize_name(name): seq.upper() for name, seq, _ in mappy.fastx_read(str(path_fasta))}


#################################################
# Export fasta files as single-FASTA format
#################################################


def export_fasta_files(TEMPDIR: Path, FASTA_ALLELES: dict, NAME: str) -> None:
    """Save multiple FASTAs in separate single-FASTA format files."""
    for identifier, sequence in FASTA_ALLELES.items():
        identifier = sanitize_name(identifier)
        contents = "\n".join([">" + identifier, sequence]) + "\n"
        path_output_fasta = Path(TEMPDIR, NAME, "fasta", f"{identifier}.fasta")
        path_output_fasta.write_text(contents)


#################################################
# save_concatenated_fastx
#################################################


def extract_extention(path_file: str | Path) -> str:
    suffixes = Path(path_file).suffixes
    return "".join(suffixes)


def is_gzip_file(path_file: str | Path) -> bool:
    """Check if a file is a GZip compressed file."""
    try:
        with Path(path_file).open("rb") as f:
            return f.read(2) == b"\x1f\x8b"
    except OSError:
        return False


def save_bam_files_to_single_fastq(TEMPDIR: Path, path_bam_files: list[Path], sample_name: str) -> None:
    path_bam_files = [str(path) for path in path_bam_files]
    path_concatenated_bam = Path(TEMPDIR, sample_name, "fastq", f"tmp_{sample_name}_{str(uuid.uuid4())}.bam")
    pysam.merge("-f", "-o", str(path_concatenated_bam), *path_bam_files)

    path_concatenated_fastq = Path(TEMPDIR, sample_name, "fastq", f"{sample_name}.fastq.gz")
    _ = pysam.fastq("-0", str(path_concatenated_fastq), str(path_concatenated_bam))

    path_concatenated_bam.unlink()


def save_fastx_files_to_single_fastq(TEMPDIR: Path, path_input_files: list[Path], sample_name: str) -> None:
    """Merge gzip and non-gzip files into a single gzip file."""
    with gzip.open(Path(TEMPDIR, sample_name, "fastq", f"{sample_name}.fastq.gz"), "wb") as merged_file:
        for path_file in path_input_files:
            if is_gzip_file(path_file):
                with gzip.open(path_file, "rb") as f:
                    merged_file.write(f.read())
            else:
                with open(path_file) as f:
                    merged_file.write(f.read().encode())


def save_inputs_as_single_fastx(TEMPDIR: Path, path_directory: Path, sample_name: str) -> None:
    file_suffix = {".fa", ".fq", ".fasta", ".fastq", ".fa.gz", ".fq.gz", ".fasta.gz", ".fastq.gz", ".bam"}
    path_input_files = [path for path in path_directory.iterdir() if extract_extention(path) in file_suffix]
    if all(extract_extention(path) == ".bam" for path in path_directory.iterdir()):
        save_bam_files_to_single_fastq(TEMPDIR, path_input_files, sample_name)
    else:
        save_fastx_files_to_single_fastq(TEMPDIR, path_input_files, sample_name)


#################################################
# save_subsetted_fastx
#################################################


def read_lines(path_file) -> list[str]:
    if is_gzip_file(path_file):
        with gzip.open(path_file, "rt") as f:
            return [line.strip() for line in f]
    else:
        with open(path_file) as f:
            return [line.strip() for line in f]


def parse_fastq(fastq) -> list[dict]:
    fastq_parsed = []
    iterator = iter(fastq)
    try:
        while True:
            header = next(iterator).strip()
            sequence = next(iterator).strip()
            annotate = next(iterator).strip()
            quality = next(iterator).strip()
            fastq_parsed.append({"header": header, "sequence": sequence, "annotate": annotate, "quality": quality})
    except StopIteration:
        pass
    return fastq_parsed


def save_subsetted_fastx(path_fastq: str | Path, num_reads: int = 10_000) -> None:
    """If the number of control reads is too high, it unnecessarily slows down the computation speed. Therefore, we perform random sampling to reduce the number of reads to below 10,000."""

    fastq: list[str] = read_lines(path_fastq)
    if len(fastq) // 4 <= 10_000:
        return None

    fastq: list[dict] = parse_fastq(fastq)

    random.seed(1)
    fastq_subset = random.sample(fastq, num_reads)
    with gzip.open(path_fastq, "wt") as f:
        for read in fastq_subset:
            f.write(f"{read['header']}\n{read['sequence']}\n{read['annotate']}\n{read['quality']}\n")
