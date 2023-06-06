from __future__ import annotations

import re
from itertools import groupby


def split_cigar(CIGAR: str) -> list[str]:
    cigar = re.split(r"([MDISH=X])", CIGAR)
    n = len(cigar)
    cigar_split = []
    for i, j in zip(range(0, n, 2), range(1, n, 2)):
        cigar_split.append(cigar[i] + cigar[j])
    return cigar_split


def trim_softclip(CIGAR: str, SEQ: str) -> str:
    cigar_split = split_cigar(CIGAR)
    if cigar_split[0].endswith("S"):
        SEQ = SEQ[int(cigar_split[0][:-1]) :]
    if cigar_split[-1].endswith("S"):
        SEQ = SEQ[: -int(cigar_split[-1][:-1])]
    return SEQ


def split_contents(sam):
    sam_headers = [s for s in sam if s[0].startswith("@")]
    sam_contents = [s for s in sam if not s[0].startswith("@") and s[9] != "*"]
    sam_contents.sort(key=lambda x: [x[0], int(x[3])])
    return sam_headers, sam_contents


def check_microhomology(current_seq_trimmed, next_seq_trimmed):
    len_microhomology = 0
    for i in range(1, min(len(current_seq_trimmed), len(next_seq_trimmed))):
        if current_seq_trimmed[-i:] == next_seq_trimmed[:i][::-1]:
            len_microhomology = i
    return len_microhomology


def format_next_align(next_align, len_microhomology, next_cigar, next_seq_trimmed, next_qual_trimmed):
    next_align[3] = str(int(next_align[3]) + len_microhomology)
    next_cigar_split = [c for c in split_cigar(next_cigar) if not re.search(r"[SH]$", c)]
    next_cigar_split[0] = str(int(next_cigar_split[0][:-1]) - len_microhomology) + next_cigar_split[0][-1]
    # TODO Check the condition when the beggining of CIGAR < 0
    if "-" in next_cigar_split[0]:
        return None
    next_align[5] = "".join(next_cigar_split)
    next_align[9] = next_seq_trimmed[len_microhomology:]
    next_align[10] = next_qual_trimmed[len_microhomology:]
    return next_align


def process_alignments(alignments):
    idx = 0
    while idx < len(alignments) - 1:
        current_align, next_align = alignments[idx], alignments[idx + 1]
        #
        current_cigar, next_cigar = current_align[5], next_align[5]
        current_seq, next_seq = current_align[9], next_align[9]
        current_qual, next_qual = current_align[10], next_align[10]
        # trim softclip
        current_seq_trimmed = trim_softclip(current_cigar, current_seq)
        next_seq_trimmed = trim_softclip(next_cigar, next_seq)
        next_qual_trimmed = trim_softclip(next_cigar, next_qual)
        # Check the length of microhomology
        len_microhomology = check_microhomology(current_seq_trimmed, next_seq_trimmed)
        # Continue if no microhomology exists
        if len_microhomology == 0:
            idx += 1
            continue
        # Correct start, cigar, seq, qual if microhomology exists
        next_align = format_next_align(next_align, len_microhomology, next_cigar, next_seq_trimmed, next_qual_trimmed)
        if next_align is not None:
            alignments[idx + 1] = next_align
        idx += 1


def remove_microhomology(sam):
    sam_headers, sam_contents = split_contents(sam)
    sam_trimmed = sam_headers.copy()
    for _, group in groupby(sam_contents, key=lambda x: x[0]):
        alignments = list(group)
        if len(alignments) == 1:
            sam_trimmed.append(alignments[0])
            continue
        process_alignments(alignments)
        sam_trimmed.extend(alignments)
    return sam_trimmed
