from __future__ import annotations

import gzip
import random
import re
import uuid
from collections.abc import Iterator
from itertools import islice
from pathlib import Path

import pysam

from DAJIN2.utils.fileio import detect_fastx_format, is_gzip_file, read_fasta, read_fastq, sanitize_name, write_fastq

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
    return {sanitize_name(record["identifier"]): record["sequence"].upper() for record in read_fasta(path_fasta)}


#################################################
# Export fasta files as single-FASTA format
#################################################


def export_fasta_files(ARGS, is_control: bool = False) -> None:
    """Save multiple FASTAs in separate single-FASTA format files."""
    tempdir = ARGS.tempdir
    if is_control:
        sample_name = ARGS.control_name
    else:
        sample_name = ARGS.sample_name

    for identifier, sequence in ARGS.fasta_alleles.items():
        identifier = sanitize_name(identifier)
        contents = "\n".join([">" + identifier, sequence]) + "\n"
        path_output_fasta = Path(tempdir, sample_name, "fasta", f"{identifier}.fasta")
        path_output_fasta.write_text(contents)


#################################################
# save_concatenated_fastx
#################################################


def extract_extention(path_file: str | Path) -> str:
    suffixes = Path(path_file).suffixes
    return "".join(suffixes)


def convert_bam_to_fastq(path_bam: str | Path, path_fastq: str | Path, quality_score: str = "I"):
    """Convert a BAM file to a gzipped FASTQ file."""
    with (
        pysam.AlignmentFile(str(path_bam), "rb", check_sq=False) as bam_file,
        gzip.open(path_fastq, "wt") as fastq_file,
    ):
        for read in bam_file:
            if read.qual is None:  # FASTA format
                read.qual = quality_score * len(read.query_sequence)
            fastq_entry = f"@{read.query_name}\n{read.query_sequence}\n+\n{read.qual}\n"
            fastq_file.write(fastq_entry)


def save_bam_files_to_single_fastq(TEMPDIR: Path, path_bam_files: list[Path], sample_name: str) -> None:
    path_bam_files = [str(path) for path in path_bam_files]
    path_concatenated_bam = Path(TEMPDIR, sample_name, "fastq", f"tmp_{sample_name}_{str(uuid.uuid4())}.bam")
    pysam.merge("-f", "-o", str(path_concatenated_bam), *path_bam_files)

    path_concatenated_fastq = Path(TEMPDIR, sample_name, "fastq", f"{sample_name}.fastq.gz")
    convert_bam_to_fastq(str(path_concatenated_bam), path_concatenated_fastq)

    path_concatenated_bam.unlink()


def merge_fastx_files_to_single_fastx(TEMPDIR: Path, path_input_files: list[Path], sample_name: str) -> None:
    """Merge gzip and non-gzip FASTX files into a single gzip FASTX file."""
    path_concatenated_fastx = Path(TEMPDIR, sample_name, "fastq", f"{sample_name}.fastx.gz")
    with gzip.open(path_concatenated_fastx, "wt") as output_file:
        for path_input_file in path_input_files:
            if is_gzip_file(path_input_file):
                with gzip.open(path_input_file, "rt") as input_file:
                    for line in input_file:
                        output_file.write(line)
            else:
                with open(path_input_file) as input_file:
                    for line in input_file:
                        output_file.write(line)


def convert_fasta_to_fastq(path_fasta: str | Path, path_fastq: str | Path, quality_score: str = "I"):
    """Convert a FASTA file to a gzipped FASTQ file with a default quality score."""
    fasta: Iterator[dict[str, str]] = read_fasta(path_fasta)
    with gzip.open(path_fastq, "wt") as fastq_file:
        for record in fasta:
            identifier = f"@{record['identifier']}"
            sequence = record["sequence"]
            separator = "+"
            quality = quality_score * len(sequence)
            fastq_file.write(f"{identifier}\n{sequence}\n{separator}\n{quality}\n")


def save_fastx_files_to_single_fastq(TEMPDIR: Path, path_input_files: list[Path], sample_name: str) -> None:
    merge_fastx_files_to_single_fastx(TEMPDIR, path_input_files, sample_name)
    path_concatenated_fastx = Path(TEMPDIR, sample_name, "fastq", f"{sample_name}.fastx.gz")
    path_concatenated_fastq = Path(TEMPDIR, sample_name, "fastq", f"{sample_name}.fastq.gz")
    if detect_fastx_format == "FASTA":
        convert_fasta_to_fastq(path_concatenated_fastx, path_concatenated_fastq)
        path_concatenated_fastx.unlink()
    else:
        path_concatenated_fastx.rename(path_concatenated_fastq)
    path_concatenated_fastx.unlink(missing_ok=True)


def _get_path_input_files(path_directory: Path, file_suffix: set[str]) -> list[Path]:
    path_input_files = []
    for path in path_directory.iterdir():
        ext = extract_extention(path)
        if ext.endswith(".fai") or ext.endswith(".bai"):
            continue
        if ext in file_suffix:
            path_input_files.append(path)
    return path_input_files


def save_inputs_as_single_fastq(ARGS, is_control: bool = False) -> None:
    tempdir = ARGS.tempdir
    if is_control:
        path_directory = ARGS.path_control
        sample_name = ARGS.control_name
    else:
        path_directory = ARGS.path_sample
        sample_name = ARGS.sample_name

    file_suffix = {".fa", ".fq", ".fasta", ".fastq", ".fa.gz", ".fq.gz", ".fasta.gz", ".fastq.gz", ".bam"}
    path_input_files = _get_path_input_files(path_directory, file_suffix)

    if path_input_files and all(extract_extention(path) == ".bam" for path in path_input_files):
        save_bam_files_to_single_fastq(tempdir, path_input_files, sample_name)
    else:
        save_fastx_files_to_single_fastq(tempdir, path_input_files, sample_name)


#################################################
# overwrite_with_downsampled_fastq
#################################################


def is_iterator_length_below_limits(iterator: Iterator, num_limits: int):
    return sum(1 for _ in islice(iterator, num_limits + 1)) <= num_limits


def overwrite_with_downsampled_fastq(path_fastq: str | Path, num_reads=10_000) -> None:
    """If the number of control reads is too high, it unnecessarily slows down the computation speed. Therefore, we perform random sampling to reduce the number of reads to below 10,000."""

    random.seed(1)

    num_limits = num_reads * 4

    if is_iterator_length_below_limits(read_fastq(path_fastq), num_limits):
        return None

    reads = list(read_fastq(path_fastq))
    sampled_reads = random.sample(reads, num_reads)
    write_fastq(sampled_reads, path_fastq, is_gzip=True)
