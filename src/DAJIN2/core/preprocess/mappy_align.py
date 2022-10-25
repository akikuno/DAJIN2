from __future__ import annotations
from collections.abc import Generator
from copy import deepcopy
from pathlib import Path
from typing import Union
import cstag
import mappy


def revcomp(sequence: str) -> str:
    complement = {"A": "T", "C": "G", "G": "C", "T": "A"}
    return "".join(complement[nt] for nt in sequence[::-1])


def to_sam(path_reference_fasta: str, path_query_fastq: str, cslong: bool = True) -> Generator[str]:
    """Align seqences using mappy and Convert PAF to SAM

    Args:
        path_reference_fasta (str): Path of reference fasta
        path_query_fastq (str): Path of query fasta/fastq
        cslong (bool, optional): long formatted CS tag if True. Defaults to True

    Returns:
        list: List of SAM
    """
    # SQ header
    SAM = [f"@SQ\tSN:{n}\tLN:{len(s)}" for n, s, _ in mappy.fastx_read(path_reference_fasta)]
    # Mappy
    ref = mappy.Aligner(path_reference_fasta)
    if not ref:
        raise AttributeError(f"Failed to load {path_reference_fasta}")
    for MAPPY_NAME, MAPPY_SEQ, MAPPY_QUAL in mappy.fastx_read(path_query_fastq):
        for hit in ref.map(MAPPY_SEQ, cs=True):
            query_seq = deepcopy(MAPPY_SEQ)
            query_qual = deepcopy(MAPPY_QUAL)
            # flag
            if hit.is_primary:
                flag = 0 if hit.strand == 1 else 16
            else:
                flag = 2048 if hit.strand == 1 else 2064
            # Append softclips to CIGAR
            cigar = hit.cigar_str
            if hit.q_st > 0:
                softclip = str(hit.q_st) + "S"
                cigar = softclip + cigar if hit.strand == 1 else cigar + softclip
            if len(MAPPY_SEQ) - hit.q_en > 0:
                softclip = str(len(MAPPY_SEQ) - hit.q_en) + "S"
                cigar = cigar + softclip if hit.strand == 1 else softclip + cigar
            # Revcomp
            if hit.strand == -1:
                query_seq = revcomp(query_seq)
                query_qual = query_qual[::-1]
            query_seq = query_seq.upper()
            # cslong
            cs = "cs:Z:" + hit.cs
            if cslong:
                cs = cstag.lengthen(hit.cs, cigar, query_seq)
            # summarize
            alignment = [
                MAPPY_NAME,
                flag,
                hit.ctg,
                hit.r_st + 1,
                hit.mapq,
                cigar,
                "*",
                0,
                0,
                query_seq,
                query_qual,
                cs,
            ]
            alignment = [str(a) for a in alignment]
            SAM.append("\t".join(alignment))
    for record in SAM:
        yield record


def output_sam(TEMPDIR, path_fasta, name_fasta, path_fastq, name_fastq):
    sam = to_sam(str(path_fasta), path_fastq)
    output_sam = Path(TEMPDIR, "sam", f"{name_fastq}_{name_fasta}.sam")
    output_sam.write_text("\n".join(sam))


########################################################################
# Create faidx
########################################################################


def make_faidx(path_fasta: Union[Path, str]) -> str:
    fasta = Path(path_fasta).read_text().split()
    name = fasta[0].strip(">")
    length = len("".join(fasta[1:]))
    offset = len(fasta[0]) + 1
    linebase = len(fasta[1])
    linewidth = len(fasta[1]) + 1
    faidx = list(map(str, [name, length, offset, linebase, linewidth]))
    faidx = "\t".join(faidx) + "\n"
    return faidx
