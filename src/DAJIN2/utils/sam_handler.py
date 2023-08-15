from __future__ import annotations

import re


def split_cigar(CIGAR: str) -> list[str]:
    cigar = re.split(r"([MIDNSH=X])", CIGAR)
    return [cigar[i] + cigar[j] for i, j in zip(range(0, len(cigar), 2), range(1, len(cigar), 2))]


def calculate_alignment_length(CIGAR: str) -> int:
    alignment_length = 0
    for c in split_cigar(CIGAR):
        if re.search(r"[MDN=X]", c[-1]):
            alignment_length += int(c[:-1])
    return alignment_length
