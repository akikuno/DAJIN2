from __future__ import annotations
import re
from collections import defaultdict


def revcomp(sequence: str) -> str:
    complement = {"A": "T", "C": "G", "G": "C", "T": "A"}
    return "".join(complement[nt] for nt in sequence[::-1])


def split_cigar(CIGAR: str) -> list[str]:
    cigar = re.split(r"([MDISH=X])", CIGAR)
    n = len(cigar)
    cigar_split = []
    for i, j in zip(range(0, n, 2), range(1, n, 2)):
        cigar_split.append(cigar[i] + cigar[j])
    return cigar_split


def calc_length(CIGAR: str) -> int:
    cigar = split_cigar(CIGAR)
    seq_length = 0
    for c in cigar:
        if re.search(r"[MD=X]", c[-1]):
            seq_length += int(c[:-1])
    return seq_length


sam_flags = [str(s) for s in [0, 16, 2048, 2064]]


def reverse_sam(sam_contents: list[str], genome_end: int) -> list[str]:
    sam_reversed = []
    for sam_content in sam_contents:
        sam_flag = sam_content[1]
        if sam_flag == sam_flags[0]:
            sam_flag = sam_flags[1]
        elif sam_flag == sam_flags[1]:
            sam_flag = sam_flags[0]
        elif sam_flag == sam_flags[2]:
            sam_flag = sam_flags[3]
        else:
            sam_flag = sam_flags[2]
        sam_content[1] = sam_flag
        sam_cigar = sam_content[5]
        sam_cigar = "".join(split_cigar(sam_cigar)[::-1])
        sam_content[5] = sam_cigar
        sam_start = int(sam_content[3])
        sam_length = calc_length(sam_cigar)
        sam_content[3] = str(genome_end - (sam_start + sam_length) + 2)
        sam_content[9] = revcomp(sam_content[9])
        sam_content[10] = sam_content[10][::-1]
        sam_reversed.append(sam_content)
    return sam_reversed


def realign(sam_header: list[str], sam_contents: list[str], genome_coodinates: dict, chrom_size: int) -> list[str]:
    for s in sam_header:
        if s[0] != "@SQ":
            continue
        s[1] = f'SN:{genome_coodinates["chr"]}'
        s[2] = f"LN:{chrom_size}"
    for s in sam_contents:
        s[2] = genome_coodinates["chr"]
    if genome_coodinates["strand"] == "-":
        sam_contents = reverse_sam(sam_contents, genome_coodinates["end"])
    else:
        s[3] = str(int(s[3]) + genome_coodinates["start"] - 1)
    sam_realined = sam_header + sam_contents
    sam_realined = ["\t".join(s) for s in sam_realined]
    return "\n".join(sam_realined)


def group_by_name(sam_contents: list[str], clust_sample: list[dict]) -> dict[list]:
    sam_contents.sort()
    clust_sample_qname = sorted(clust_sample, key=lambda x: x["QNAME"])
    clust_sample_qname_set = set()
    for qnames in clust_sample_qname:
        qname = qnames["QNAME"]
        clust_sample_qname_set.add(qname)
    sam_groups = defaultdict(list)
    idx_left = 0
    idx_right = 0
    while idx_left < len(sam_contents) and idx_right < len(clust_sample_qname):
        read_left = sam_contents[idx_left][:-1]
        read_right = clust_sample_qname[idx_right]
        qname_left = read_left[0]
        qname_right = read_right["QNAME"]
        if qname_left not in clust_sample_qname_set:
            idx_left += 1
            continue
        if qname_left == qname_right:
            key = read_right["NAME"]
            sam_groups[key].append(read_left)
            idx_left += 1
        else:
            idx_right += 1
    return sam_groups

