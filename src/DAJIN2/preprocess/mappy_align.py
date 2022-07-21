from __future__ import annotations
from collections.abc import Generator

import re
from itertools import groupby
from collections import deque

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
        cslong (bool, optional): long formatted CS tag if True. Defaults to True.

    Returns:
        list: List of SAM
    """
    # SQ header
    SAM = [f"@SQ\tSN:{n}\tLN:{len(s)}" for n, s, _ in mappy.fastx_read(path_reference_fasta)]
    # Mappy
    ref = mappy.Aligner(path_reference_fasta)
    if not ref:
        raise AttributeError(f"Failed to load f{path_reference_fasta}")
    for query_name, query_sequence, query_quality in mappy.fastx_read(path_query_fastq):
        for hit in ref.map(query_sequence, cs=True):
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
            if len(query_sequence) - hit.q_en > 0:
                softclip = str(len(query_sequence) - hit.q_en) + "S"
                cigar = cigar + softclip if hit.strand == 1 else softclip + cigar
            # Revcomp
            if not hit.strand == 1:
                query_sequence = revcomp(query_sequence)
            query_sequence = query_sequence.upper()
            # cslong
            cs = "cs:Z:" + hit.cs
            if cslong:
                cs = cstag.lengthen(hit.cs, cigar, query_sequence)
            # summarize
            alignment = [
                query_name,
                flag,
                hit.ctg,
                hit.r_st + 1,
                hit.mapq,
                cigar,
                "*",
                0,
                0,
                query_sequence,
                query_quality,
                cs,
            ]
            alignment = [str(a) for a in alignment]
            SAM.append("\t".join(alignment))
    for record in SAM:
        yield record


# def remove_unmapped_reads(sam: list[str]) -> Generator[str]:
#     sam_mapped_reads = []
#     for record in sam:
#         if record.startswith("@"):
#             sam_mapped_reads.append(record)
#             continue
#         if not record.split("\t")[2] == "*":
#             sam_mapped_reads.append(record)
#     for record in sam_mapped_reads:
#         yield record


# def remove_overlapped_reads(sam: list[str]) -> Generator[str]:
#     sam = [s.split("\t") for s in sam]
#     sam.sort(key=lambda x: x[0])
#     sam_groupby = groupby(sam, lambda x: x[0])
#     sam_nonoverlapped = deque()
#     for key, record_gropby in sam_groupby:
#         if key.startswith("@"):
#             for record in record_gropby:
#                 sam_nonoverlapped.appendleft("\t".join(record))
#             continue
#         records = sorted(record_gropby, key=lambda x: x[3])
#         is_overraped = False
#         end_of_previous_read = -1
#         for record in records:
#             start_of_current_read = int(record[3])
#             if end_of_previous_read > start_of_current_read:
#                 is_overraped = True
#                 break
#             record_length = 0
#             cigar = record[5]
#             cigar_split = re.split(r"([A-Z])", cigar)
#             for i, cigar in enumerate(cigar_split):
#                 if cigar == "M" or cigar == "D":
#                     record_length += int(cigar_split[i - 1])
#             end_of_previous_read = start_of_current_read + record_length - 1
#         if is_overraped:
#             continue
#         for record in records:
#             sam_nonoverlapped.append("\t".join(record))
#     for record in list(sam_nonoverlapped):
#         yield record

