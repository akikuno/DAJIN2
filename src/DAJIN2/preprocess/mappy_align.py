from __future__ import annotations

import re
from itertools import groupby
import cstag
import mappy


def revcomp(sequence: str) -> str:
    complement = {"A": "T", "C": "G", "G": "C", "T": "A"}
    return "".join(complement[nt] for nt in sequence[::-1])


def to_sam(path_reference_fasta: str, path_query_fastq: str, cslong: bool = True) -> list[str]:
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
    return SAM


def remove_unmapped_reads(sam: list[str]) -> list[str]:
    sam_mapped_reads = []
    for record in sam:
        if record.startswith("@"):
            sam_mapped_reads.append(record)
            continue
        if not record.split("\t")[2] == "*":
            sam_mapped_reads.append(record)
    return sam_mapped_reads


def remove_overlapped_reads(sam: list[str]) -> list[str]:
    sam = [s.split("\t") for s in sam]
    sam.sort(key=lambda x: x[0])
    sam_groupby = groupby(sam, lambda x: x[0])
    sam_nonoverlapped = []
    for key, record_gropby in sam_groupby:
        if key.startswith("@"):
            for record in record_gropby:
                sam_nonoverlapped.append("\t".join(record))
            continue
        records = sorted(record_gropby, key=lambda x: x[3])
        is_overraped = False
        end_of_previous_read = -1
        for record in records:
            start_of_current_read = int(record[3])
            if end_of_previous_read > start_of_current_read:
                is_overraped = True
                break
            record_length = 0
            cigar = record[5]
            cigar_split = re.split(r"([A-Z])", cigar)
            for i, cigar in enumerate(cigar_split):
                if cigar == "M" or cigar == "I":
                    record_length += int(cigar_split[i - 1])
            end_of_previous_read = start_of_current_read + record_length
        if is_overraped:
            continue
        for record in records:
            sam_nonoverlapped.append("\t".join(record))
    return sam_nonoverlapped


# def remove_long_softclipped_reads(sam: list[str]) -> list[str]:
#     """Remove reads with soft clips longer than 1/10 of the reference sequence length
#     """
#     sam_short_softcliped_reads = []
#     sqheaders = dict()
#     for record in sam:
#         if record.startswith("@"):
#             sam_short_softcliped_reads.append(record)
#         if record.startswith("@SQ"):
#             for sqheader in record.split("\t"):
#                 if sqheader.startswith("SN:"):
#                     SN = sqheader.replace("SN:", "")
#                 if sqheader.startswith("LN:"):
#                     LN = int(sqheader.replace("LN:", ""))
#             sqheaders.update({SN: LN})
#             continue
#         rname = record.split("\t")[2]
#         reference_sequence_length = sqheaders[rname]
#         cigar = record.split("\t")[5]
#         cigar_split = re.split(r"([A-Z])", cigar)
#         softclip = 0
#         for i, s in enumerate(cigar_split):
#             if s == "S":
#                 softclip += int(cigar_split[i - 1])
#         if softclip < reference_sequence_length // 10:
#             sam_short_softcliped_reads.append(record)
#     return sam_short_softcliped_reads

