from __future__ import annotations

import cstag
import mappy

from pathlib import Path
from typing import Generator

from DAJIN2.utils.dna_handler import revcomp


def to_sam(
    path_reference_fasta: Path, path_query_fastx: Path, preset: str = "map-ont", threads: int = 1, cslong: bool = True
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
    path_reference_fasta = str(path_reference_fasta)
    path_query_fastx = str(path_query_fastx)

    SAM = [f"@SQ\tSN:{n}\tLN:{len(s)}" for n, s, _ in mappy.fastx_read(path_reference_fasta)]

    ref = mappy.Aligner(path_reference_fasta, preset=preset, n_threads=threads)
    if not ref:
        raise ValueError(f"Failed to load {path_reference_fasta}")

    for QUERY_NAME, QUERY_SEQ, QUERY_QUAL in mappy.fastx_read(path_query_fastx):
        for hit in ref.map(QUERY_SEQ, cs=True):
            query_seq = QUERY_SEQ.upper()
            query_qual = QUERY_QUAL

            # Report flag
            if hit.is_primary:
                flag = 0 if hit.strand == 1 else 16
            else:
                flag = 2048 if hit.strand == 1 else 2064

            # Handle reverse complement for negative strand
            if hit.strand == -1:
                query_seq = revcomp(query_seq)
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

    for record in SAM:
        yield record


########################################################################
# main
########################################################################


def generate_sam(ARGS, paths_fasta: list[str], is_control: bool = False, is_insertion: bool = False) -> None:
    if is_control:
        path_fastq = ARGS.path_control
        name = ARGS.control_name
    else:
        path_fastq = ARGS.path_sample
        name = ARGS.sample_name

    out_directory = Path(ARGS.tempdir, name, "sam")

    for path_fasta in paths_fasta:
        name_fasta = Path(path_fasta).stem
        for preset in ["map-ont", "splice"]:
            sam = to_sam(path_fasta, path_fastq, preset=preset, threads=ARGS.threads)
            if is_control and is_insertion:
                path_sam = Path(out_directory, f"{preset}_{name_fasta}_{ARGS.sample_name}.sam")
            else:
                path_sam = Path(out_directory, f"{preset}_{name_fasta}.sam")
            path_sam.write_text("\n".join(sam))


########################################################################
# Create faidx
########################################################################


def make_faidx(path_fasta: str | Path) -> str:
    fasta = Path(path_fasta).read_text().split()
    name, length, offset = fasta[0].strip(">"), len("".join(fasta[1:])), len(fasta[0]) + 1
    linebase, linewidth = len(fasta[1]), len(fasta[1]) + 1
    return "\t".join(map(str, [name, length, offset, linebase, linewidth])) + "\n"
