from __future__ import annotations

import re


def revcomp(sequence: str) -> str:
    """Reverse complement a DNA sequence."""
    complement = {"A": "T", "C": "G", "G": "C", "T": "A"}
    return "".join(complement[nt] for nt in sequence[::-1])


def split_cigar(CIGAR: str) -> list[str]:
    """Split a CIGAR string into its individual elements."""
    cigar = re.split(r"([MDISH=X])", CIGAR)
    n = len(cigar)
    cigar_split = []
    for i, j in zip(range(0, n, 2), range(1, n, 2)):
        cigar_split.append(cigar[i] + cigar[j])
    return cigar_split


def calc_length(CIGAR: str) -> int:
    """Calculate the length of the sequence represented by a CIGAR string."""
    cigar = split_cigar(CIGAR)
    return sum(int(c[:-1]) for c in cigar if c[-1] in "MD=X")


def reverse_sam(sam_contents: list[list[str]], genome_end: int) -> list[str]:
    """Reverse and complement SAM entries."""
    flag_map = {"0": "16", "16": "0", "2048": "2064", "2064": "2048"}
    sam_reversed = []
    for sam_content in sam_contents:
        sam_update = sam_content.copy()
        sam_flag = sam_content[1]
        sam_update[1] = flag_map.get(sam_flag, sam_flag)
        sam_cigar = "".join(split_cigar(sam_content[5])[::-1])
        sam_update[5] = sam_cigar
        sam_start = int(sam_content[3])
        sam_length = calc_length(sam_cigar)
        sam_update[3] = str(genome_end - (sam_start + sam_length) + 2)
        sam_update[9] = revcomp(sam_content[9])
        sam_update[10] = sam_content[10][::-1]
        sam_reversed.append(sam_update)
    return sam_reversed
