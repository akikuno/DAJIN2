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


def split_contents(sam: list[list[str]]) -> list[list[str]]:
    sam_headers = [s for s in sam if s[0].startswith("@")]
    sam_contents = [s for s in sam if not s[0].startswith("@") and s[9] != "*"]
    sam_contents.sort(key=lambda x: [x[0], int(x[3])])
    return sam_headers, sam_contents


def check_microhomology(curr_seq_trimmed: str, next_seq_trimmed: str) -> int:
    len_microhomology = 0
    for i in range(1, min(len(curr_seq_trimmed), len(next_seq_trimmed))):
        if curr_seq_trimmed[-i:] == next_seq_trimmed[:i]:
            len_microhomology = i
    return len_microhomology


def count_mutations_in_microhomology(cigar_split: list, len_microhomology: int) -> int:
    """Count the number of mutations within the microhomology based on a CIGAR list."""
    total_num = 0
    count_mutation = 0
    for cigar in cigar_split:
        num, op = int(cigar[:-1]), cigar[-1]
        total_num += num
        # Stop counting when total_num reaches the length of microhomology
        if total_num >= len_microhomology:
            break
        # Increment mutation count for non-match operations
        if op != "M":
            count_mutation += num
    return count_mutation


def trim_cigar_on_microhomology(cigar_split: list, len_microhomology: int) -> tuple(str, int):
    """Trim CIGAR based on the length of the microhomology."""
    cigar_trimmed = cigar_split[::-1]
    total_num = 0
    num_deleletion = 0
    while total_num < len_microhomology:
        cigar = cigar_trimmed.pop()
        num, op = int(cigar[:-1]), cigar[-1]
        # Skip deletion operations
        if op == "D":
            num_deleletion += num
            continue
        total_num += num
        # Only append to the list when total_num is greater or equal to len_microhomology
        if total_num > len_microhomology:
            num = total_num - len_microhomology
            cigar_trimmed.append(f"{num}{op}")
    return "".join(cigar_trimmed[::-1]), num_deleletion


def process_alignments(alignments: list[list[str]]) -> list[list[str]]:
    idx = 0
    while idx < len(alignments) - 1:
        curr_align, next_align = alignments[idx], alignments[idx + 1]
        curr_cigar, next_cigar = curr_align[5], next_align[5]
        curr_seq, next_seq = curr_align[9], next_align[9]
        next_qual = next_align[10]
        # trim softclip
        curr_seq_trimmed = trim_softclip(curr_cigar, curr_seq)
        next_seq_trimmed = trim_softclip(next_cigar, next_seq)
        next_qual_trimmed = trim_softclip(next_cigar, next_qual)
        # Check the length of microhomology
        len_microhomology = check_microhomology(curr_seq_trimmed, next_seq_trimmed)
        # Continue if no microhomology exists
        if len_microhomology == 0:
            idx += 1
            continue
        # Update CIGAR
        curr_cigar_splitted = [c for c in split_cigar(curr_cigar) if not re.search(r"[SH]$", c)]
        next_cigar_splitted = [c for c in split_cigar(next_cigar) if not re.search(r"[SH]$", c)]
        mutations_in_curr = count_mutations_in_microhomology(curr_cigar_splitted[::-1], len_microhomology)
        mutations_in_next = count_mutations_in_microhomology(next_cigar_splitted, len_microhomology)
        # Update START, CIGAR, SEQ, QUAL to the more mismatched alignment within microhomology
        if mutations_in_next >= mutations_in_curr:
            alignments[idx + 1][5], num_del = trim_cigar_on_microhomology(next_cigar_splitted, len_microhomology)
            alignments[idx + 1][3] = str(int(next_align[3]) + len_microhomology + num_del)
            alignments[idx + 1][9] = next_seq_trimmed[len_microhomology:]
            alignments[idx + 1][10] = next_qual_trimmed[len_microhomology:]
        else:
            alignments[idx][5], num_del = trim_cigar_on_microhomology(next_cigar_splitted, len_microhomology)
            alignments[idx][3] = str(int(alignments[idx][3]) + num_del)
            alignments[idx][9] = next_seq_trimmed[:-len_microhomology]
            alignments[idx][10] = next_qual_trimmed[:-len_microhomology:]
        idx += 1


def remove_microhomology(sam: list[list[str]]) -> list[list[str]]:
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
