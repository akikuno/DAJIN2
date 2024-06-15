from __future__ import annotations

from pathlib import Path
from typing import Generator

import cstag
import mappy

from DAJIN2.utils import dna_handler


def to_sam(
    path_reference_fasta: Path,
    path_query_fastx: Path,
    preset: str = "map-ont",
    threads: int = 1,
    options: dict = None,
    cslong: bool = True,
) -> Generator[str]:
    """Align sequences using mappy and Convert PAF to SAM.

    Args:
        path_reference_fasta (Path): Path of reference fasta.
        path_query_fastx (Path): Path of query fasta/fastq.
        preset (str, optional): Alignment preset. Defaults to "map-ont".
        threads (int, optional): Number of threads to use. Defaults to 1.
        cslong (bool, optional): Use long formatted CS tag if True. Defaults to True.

    Yields:
        str: SAM formatted alignment.
    """
    if options is None:
        options = {}

    path_reference_fasta = str(path_reference_fasta)
    path_query_fastx = str(path_query_fastx)

    SAM = [f"@SQ\tSN:{n}\tLN:{len(s)}" for n, s, _ in mappy.fastx_read(path_reference_fasta)]

    ref = mappy.Aligner(path_reference_fasta, preset=preset, n_threads=threads, **options)
    if not ref:
        raise ValueError(f"Failed to load {path_reference_fasta}")

    for QUERY_NAME, QUERY_SEQ, QUERY_QUAL in mappy.fastx_read(path_query_fastx):
        for hit in ref.map(QUERY_SEQ, cs=True):
            query_seq = QUERY_SEQ.upper()
            query_qual = QUERY_QUAL

            # Skip multi-mapping reads
            if hit.mapq == 0:
                continue

            # Report flag
            if hit.is_primary:
                flag = 0 if hit.strand == 1 else 16
            else:
                flag = 2048 if hit.strand == 1 else 2064

            # Handle reverse complement for negative strand
            if hit.strand == -1:
                query_seq = dna_handler.revcomp(query_seq)
                if query_qual:
                    query_qual = query_qual[::-1]

            # Append softclips to CIGAR
            cigar = hit.cigar_str
            if hit.q_st > 0:
                softclip = f"{hit.q_st}S"
                cigar = softclip + cigar if hit.strand == 1 else cigar + softclip
            if len(query_seq) - hit.q_en > 0:
                softclip = f"{len(query_seq) - hit.q_en}S"
                cigar = cigar + softclip if hit.strand == 1 else softclip + cigar

            # Convert to CS tag's long format
            if cslong:
                cs = cstag.lengthen(hit.cs, cigar, query_seq, prefix=True)

            # Summarize
            alignment = [
                QUERY_NAME,
                str(flag),
                hit.ctg,
                str(hit.r_st + 1),
                str(hit.mapq),
                cigar,
                "*",
                "0",
                "0",
                query_seq,
                "*" if query_qual is None else query_qual,
                cs,
            ]

            SAM.append("\t".join(alignment))

    yield from SAM


########################################################################
# main
########################################################################


def generate_sam(
    ARGS, paths_fasta: list[str], mappy_options: dict = None, is_control: bool = False, is_insertion: bool = False
) -> None:
    if mappy_options is None:
        mappy_options = {}

    if is_control:
        path_fastq = Path(ARGS.tempdir, ARGS.control_name, "fastq", f"{ARGS.control_name}.fastq.gz")
        name = ARGS.control_name
    else:
        path_fastq = Path(ARGS.tempdir, ARGS.sample_name, "fastq", f"{ARGS.sample_name}.fastq.gz")
        name = ARGS.sample_name

    for path_fasta in paths_fasta:
        name_fasta = Path(path_fasta).stem
        len_sequence = len(Path(path_fasta).read_text().split("\n")[1])
        if len_sequence < 500:
            presets = ["sr"]
        else:
            presets = ["map-ont", "splice"]

        path_sam_directory = Path(ARGS.tempdir, name, "sam", name_fasta)
        path_sam_directory.mkdir(parents=True, exist_ok=True)

        for preset in presets:
            sam = to_sam(path_fasta, path_fastq, preset=preset, threads=ARGS.threads, options=mappy_options)

            if is_control and is_insertion:
                path_sam_file = Path(path_sam_directory, f"{ARGS.sample_name}_{preset}.sam")
            else:
                path_sam_file = Path(path_sam_directory, f"{preset}.sam")

            path_sam_file.write_text("\n".join(sam))


########################################################################
# Create faidx
########################################################################


def make_faidx(path_fasta: str | Path) -> str:
    fasta = Path(path_fasta).read_text().split()
    name, length, offset = fasta[0].strip(">"), len("".join(fasta[1:])), len(fasta[0]) + 1
    linebase, linewidth = len(fasta[1]), len(fasta[1]) + 1
    return "\t".join(map(str, [name, length, offset, linebase, linewidth])) + "\n"
