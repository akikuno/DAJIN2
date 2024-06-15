from __future__ import annotations

import re
from itertools import groupby

from DAJIN2.utils.dna_handler import revcomp


def split_cigar(CIGAR: str) -> list[str]:
    cigar = re.split(r"([MIDNSH=X])", CIGAR)
    return [cigar[i] + cigar[j] for i, j in zip(range(0, len(cigar), 2), range(1, len(cigar), 2))]


def calculate_alignment_length(CIGAR: str) -> int:
    return sum(int(c[:-1]) for c in split_cigar(CIGAR) if c[-1] in "MDN=X")


def is_header(s: list[str]) -> bool:
    return s[0].startswith("@")


def is_mapped(s: list[str]) -> bool:
    return not s[0].startswith("@") and s[9] != "*"


###########################################################
# remove_overlapped_reads
###########################################################


def is_overlapping(prev_alignment: list[str], next_alignment: list[str]) -> bool:
    prev_start = int(prev_alignment[3])
    prev_alignment_len = calculate_alignment_length(prev_alignment[5])
    prev_end = prev_start + prev_alignment_len
    next_start = int(next_alignment[3])
    return prev_end >= next_start


def remove_overlapped_reads(sam: list[list[str]]) -> list[list[str]]:
    sam_headers = [s for s in sam if is_header(s)]
    sam_contents = [s for s in sam if is_mapped(s)]

    sam_contents.sort(key=lambda x: (x[0], int(x[3])))
    sam_trimmed = sam_headers.copy()

    for _, group in groupby(sam_contents, key=lambda x: x[0]):
        alignments = list(group)

        if len(alignments) == 1:
            sam_trimmed.append(alignments[0])
            continue

        if not any(is_overlapping(alignments[idx], alignments[idx + 1]) for idx in range(len(alignments) - 1)):
            sam_trimmed += alignments

    return sam_trimmed


###########################################################
# remove_microholomogy
###########################################################


def trim_softclip(cigar: str, sequence: str) -> str:
    if sequence == "*":  # In the case of FASTA format, the qual (= sequence) will be "*".
        return "*"
    cigar_split = split_cigar(cigar)
    if cigar_split[0].endswith("S"):
        sequence = sequence[int(cigar_split[0][:-1]) :]
    if cigar_split[-1].endswith("S"):
        sequence = sequence[: -int(cigar_split[-1][:-1])]
    return sequence


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


def trim_cigar_on_microhomology(cigar_split: list, len_microhomology: int) -> tuple[str, int]:
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


def trim(alignments: list[list[str]]) -> list[list[str]]:
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
            alignments[idx + 1][10] = next_qual_trimmed[len_microhomology:] if next_qual_trimmed != "*" else "*"
        else:
            alignments[idx][5], num_del = trim_cigar_on_microhomology(next_cigar_splitted, len_microhomology)
            alignments[idx][3] = str(int(alignments[idx][3]) + num_del)
            alignments[idx][9] = next_seq_trimmed[:-len_microhomology]
            alignments[idx][10] = next_qual_trimmed[:-len_microhomology:] if next_qual_trimmed != "*" else "*"
        idx += 1


def remove_microhomology(sam: list[list[str]]) -> list[list[str]]:
    sam_headers = [s for s in sam if is_header(s)]
    sam_contents = [s for s in sam if is_mapped(s)]

    sam_contents.sort(key=lambda x: (x[0], int(x[3])))
    sam_trimmed = sam_headers.copy()

    for _, group in groupby(sam_contents, key=lambda x: x[0]):
        alignments = list(group)
        if len(alignments) == 1:
            sam_trimmed.append(alignments[0])
            continue
        trim(alignments)
        sam_trimmed.extend(alignments)
    return sam_trimmed


###########################################################
# revcomp_sam
###########################################################


def reverse_flag(flag: int) -> int:
    """Reverse the flag of a SAM entry."""
    return flag | 16 if flag & 16 == 0 else flag & ~16


def revcomp_sam(sam_contents: list[list[str]], genome_end: int) -> list[str]:
    """Reverse and complement SAM entries."""
    sam_reversed = []
    for sam_content in sam_contents:
        sam_update = sam_content.copy()
        sam_update[1] = str(reverse_flag(int(sam_content[1])))
        sam_cigar = "".join(split_cigar(sam_content[5])[::-1])
        sam_update[5] = sam_cigar
        sam_start = int(sam_content[3])
        sam_length = calculate_alignment_length(sam_cigar)
        sam_update[3] = str(genome_end - (sam_start + sam_length) + 2)
        sam_update[9] = revcomp(sam_content[9])
        sam_update[10] = sam_content[10][::-1]
        sam_reversed.append(sam_update)
    return sam_reversed
