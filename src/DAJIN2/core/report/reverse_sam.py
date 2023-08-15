from __future__ import annotations

from DAJIN2.utils import sam_handler


def revcomp(sequence: str) -> str:
    """Reverse complement a DNA sequence."""
    complement = {"A": "T", "C": "G", "G": "C", "T": "A"}
    return "".join(complement[nt] for nt in sequence[::-1])


def reverse_sam(sam_contents: list[list[str]], genome_end: int) -> list[str]:
    """Reverse and complement SAM entries."""
    flag_map = {"0": "16", "16": "0", "2048": "2064", "2064": "2048"}
    sam_reversed = []
    for sam_content in sam_contents:
        sam_update = sam_content.copy()
        sam_flag = sam_content[1]
        sam_update[1] = flag_map.get(sam_flag, sam_flag)
        sam_cigar = "".join(sam_handler.split_cigar(sam_content[5])[::-1])
        sam_update[5] = sam_cigar
        sam_start = int(sam_content[3])
        sam_length = sam_handler.calculate_alignment_length(sam_cigar)
        sam_update[3] = str(genome_end - (sam_start + sam_length) + 2)
        sam_update[9] = revcomp(sam_content[9])
        sam_update[10] = sam_content[10][::-1]
        sam_reversed.append(sam_update)
    return sam_reversed
