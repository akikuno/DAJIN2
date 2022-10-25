from __future__ import annotations
import re
import random
from pathlib import Path
from typing import Union
from collections import defaultdict
from itertools import groupby
from copy import deepcopy
import pysam
import midsv


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


def write_sam(sam: list[list[str]], path_sam: Union[str, Path]):
    path_sam = str(path_sam)
    sam = ["\t".join(s) for s in sam]
    sam = "\n".join(sam)
    Path(path_sam).write_text(sam + "\n")


def remove_overlapped_reads(sam: list[list[str]]) -> list[list[str]]:
    sam_headers = [s for s in sam if s[0].startswith("@")]
    sam_contents = [s for s in sam if not s[0].startswith("@") and s[9] != "*"]
    sam_contents.sort(key=lambda x: [x[0], int(x[3])])
    sam_trimmed = sam_headers.copy()
    for _, group in groupby(sam_contents, key=lambda x: x[0]):
        alignments = list(group)
        if len(alignments) == 1:
            sam_trimmed.append(alignments[0])
            continue
        idx = 0
        flag_overlap = False
        while idx < len(alignments) - 1:
            prev_alignment = alignments[idx]
            next_alignment = alignments[idx + 1]
            prev_start = int(prev_alignment[3])
            prev_cigar = prev_alignment[5]
            prev_alignment_len = calc_length(prev_cigar)
            prev_end = prev_start + prev_alignment_len
            next_start = int(next_alignment[3])
            if prev_end >= next_start:
                flag_overlap = True
                break
            idx += 1
        if flag_overlap:
            continue
        else:
            sam_trimmed += alignments
    return sam_trimmed


def trim_softclip(CIGAR: str, SEQ: str) -> str:
    cigar_split = split_cigar(CIGAR)
    if cigar_split[0].endswith("S"):
        SEQ = SEQ[int(cigar_split[0][:-1]) :]
    if cigar_split[-1].endswith("S"):
        SEQ = SEQ[: -int(cigar_split[-1][:-1])]
    return SEQ


def remove_microhomology(sam: list[list[str]]) -> list[list[str]]:
    sam_headers = [s for s in sam if s[0].startswith("@")]
    sam_contents = [s for s in sam if not s[0].startswith("@") and s[9] != "*"]
    sam_contents.sort(key=lambda x: [x[0], int(x[3])])
    sam_trimmed = sam_headers.copy()
    for _, group in groupby(sam_contents, key=lambda x: x[0]):
        alignments = list(group)
        if len(alignments) == 1:
            sam_trimmed.append(alignments[0])
            continue
        idx = 0
        while idx < len(alignments) - 1:
            prev_align = alignments[idx]
            next_align = alignments[idx + 1]
            prev_cigar = prev_align[5]
            next_cigar = next_align[5]
            prev_seq = prev_align[9]
            next_seq = next_align[9]
            prev_qual = prev_align[10]
            next_qual = next_align[10]
            # trim softclips of sequence
            prev_seq_trimmed = trim_softclip(prev_cigar, prev_seq)
            next_seq_trimmed = trim_softclip(next_cigar, next_seq)
            # trim softclips of quality
            prev_qual_trimmed = trim_softclip(prev_cigar, prev_qual)
            next_qual_trimmed = trim_softclip(next_cigar, next_qual)
            if prev_seq_trimmed == next_seq_trimmed:
                sam_trimmed.append(prev_align)
                sam_trimmed.append(next_align)
                idx += 1
                continue
            i = 1
            len_microhomology = 0
            while i < min(len(prev_seq_trimmed), len(next_seq_trimmed)):
                if prev_seq_trimmed[-i:] == next_seq_trimmed[:i] and prev_qual_trimmed[-i:] == next_qual_trimmed[:i]:
                    len_microhomology = i
                i += 1
            if len_microhomology == 0:
                sam_trimmed.append(prev_align)
                sam_trimmed.append(next_align)
                idx += 1
                continue
            # ----------------------
            # format
            # ----------------------
            # start
            next_align[3] = str(int(next_align[3]) + len_microhomology)
            # cigar
            next_cigar_split = [c for c in split_cigar(next_cigar) if not re.search(r"[SH]$", c)]
            next_cigar_split[0] = str(int(next_cigar_split[0][:-1]) - len_microhomology) + next_cigar_split[0][-1]
            # remove reads with overlapped at softclip region
            if "-" in next_cigar_split[0]:
                idx += 1
                continue
            next_align[5] = "".join(next_cigar_split)
            # sequence
            next_align[9] = next_seq_trimmed[len_microhomology:]
            # quality
            next_align[10] = next_qual_trimmed[len_microhomology:]
            # ----------------------
            # finish
            # ----------------------
            sam_trimmed.append(prev_align)
            sam_trimmed.append(next_align)
            idx += 1
    return sam_trimmed


def reverse_sam(sam_contents: list[list[str]], genome_end: int) -> list[str]:
    sam_flags = [str(s) for s in [0, 16, 2048, 2064]]
    sam_reversed = []
    for sam_content in sam_contents:
        sam_update = sam_content.copy()
        sam_flag = sam_content[1]
        if sam_flag == sam_flags[0]:
            sam_flag = sam_flags[1]
        elif sam_flag == sam_flags[1]:
            sam_flag = sam_flags[0]
        elif sam_flag == sam_flags[2]:
            sam_flag = sam_flags[3]
        else:
            sam_flag = sam_flags[2]
        sam_update[1] = sam_flag
        sam_cigar = sam_content[5]
        sam_cigar = "".join(split_cigar(sam_cigar)[::-1])
        sam_update[5] = sam_cigar
        sam_start = int(sam_content[3])
        sam_length = calc_length(sam_cigar)
        sam_update[3] = str(genome_end - (sam_start + sam_length) + 2)
        sam_update[9] = revcomp(sam_content[9])
        sam_update[10] = sam_content[10][::-1]
        sam_reversed.append(sam_update)
    return sam_reversed


def realign(sam: list[list[str]], genome_coodinates: dict, chrom_size: int) -> list[str]:
    sam_headers = [s for s in sam if s[0].startswith("@")]
    sam_contents = [s for s in sam if not s[0].startswith("@")]
    for s in sam_headers:
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
    return sam_headers + sam_contents


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


###############################################################################
# igvjs
###############################################################################


def subset_qnames(RESULT_SAMPLE, readnum: int = 100) -> dict[set[str]]:
    qnames_by_name = defaultdict(set)
    for name, group in groupby(RESULT_SAMPLE, key=lambda x: x["NAME"]):
        group = list(group)
        qnames = [res["QNAME"] for res in group[:readnum]]
        qnames_by_name[name] = set(qnames)
    return qnames_by_name


def subset_reads(name, sam_content, qnames_by_name):
    qnames = qnames_by_name[name]
    sam_subset = [sam for sam in sam_content if sam[0] in qnames]
    return sam_subset


###############################################################################
# Output
###############################################################################


def output_bam_control(TEMPDIR, CONTROL_NAME, GENOME, GENOME_COODINATES, CHROME_SIZE, THREADS):
    randomnum = random.randint(100_000, 999_999)

    path_sam = Path(TEMPDIR, "sam", f"{CONTROL_NAME}_control.sam")
    sam = midsv.read_sam(path_sam)

    sam_update = deepcopy(sam)
    sam_update = remove_overlapped_reads(sam_update)
    sam_update = remove_microhomology(sam_update)

    if GENOME:
        sam_update = realign(sam_update, GENOME_COODINATES, CHROME_SIZE)

    path_sam_update = Path(TEMPDIR, "report", "bam", f"tmp{randomnum}_{CONTROL_NAME}_control.sam")
    write_sam(sam_update, path_sam_update)

    path_bam = Path(TEMPDIR, "report", "bam", f"{CONTROL_NAME}.bam")
    pysam.sort("-@", f"{THREADS}", "-o", str(path_bam), str(path_sam_update))
    pysam.index("-@", f"{THREADS}", str(path_bam))

    # igvjs
    sam_headers = [s for s in sam_update if s[0].startswith("@")]
    sam_contents = [s for s in sam_update if not s[0].startswith("@")]
    qnames = [s[0] for s in sam_contents[:10000]]
    qnames = set(list(set(qnames))[:100])
    sam_subset = [s for s in sam_update if s[0] in qnames]
    path_bam = Path(TEMPDIR, "report", ".igvjs", f"{CONTROL_NAME}.bam")
    write_sam(sam_headers + sam_subset, path_sam)
    pysam.sort("-@", f"{THREADS}", "-o", str(path_bam), str(path_sam))
    pysam.index("-@", f"{THREADS}", str(path_bam))

    # Remove temporary files
    sam_temp = Path(TEMPDIR, "report", "bam").glob(f"tmp{randomnum}*.sam")
    [s.unlink() for s in sam_temp]


def output_bam_sample(TEMPDIR, RESULT_SAMPLE, SAMPLE_NAME, GENOME, GENOME_COODINATES, CHROME_SIZE, THREADS):
    randomnum = random.randint(100_000, 999_999)
    path_sam = Path(TEMPDIR, "sam", f"{SAMPLE_NAME}_control.sam")
    sam = midsv.read_sam(path_sam)

    sam_update = deepcopy(sam)
    sam_update = remove_overlapped_reads(sam_update)
    sam_update = remove_microhomology(sam_update)

    if GENOME:
        sam_update = realign(sam_update, GENOME_COODINATES, CHROME_SIZE)

    path_sam_update = Path(TEMPDIR, "report", "bam", f"tmp{randomnum}_{SAMPLE_NAME}_control.sam")
    write_sam(sam_update, path_sam_update)

    path_bam = Path(TEMPDIR, "report", "bam", f"{SAMPLE_NAME}.bam")
    pysam.sort("-@", f"{THREADS}", "-o", str(path_bam), str(path_sam_update))
    pysam.index("-@", f"{THREADS}", str(path_bam))

    sam_headers = [s for s in sam_update if s[0].startswith("@")]
    sam_contents = [s for s in sam_update if not s[0].startswith("@")]
    sam_groups = group_by_name(sam_contents, RESULT_SAMPLE)
    qnames_by_name = subset_qnames(RESULT_SAMPLE)

    for name, sam_content in sam_groups.items():
        # BAM
        path_sam = Path(TEMPDIR, "report", "bam", f"tmp{randomnum}_{name}.sam")
        path_bam = Path(TEMPDIR, "report", "bam", f"{SAMPLE_NAME}_{name}.bam")
        write_sam(sam_headers + sam_content, path_sam)
        pysam.sort("-@", f"{THREADS}", "-o", str(path_bam), str(path_sam))
        pysam.index("-@", f"{THREADS}", str(path_bam))
        # igvjs
        sam_subset = subset_reads(name, sam_content, qnames_by_name)
        path_sam = Path(TEMPDIR, "report", "bam", f"tmp{randomnum}_{name}_subset.sam")
        path_bam = Path(TEMPDIR, "report", ".igvjs", f"{SAMPLE_NAME}_{name}.bam")
        write_sam(sam_headers + sam_subset, path_sam)
        pysam.sort("-@", f"{THREADS}", "-o", str(path_bam), str(path_sam))
        pysam.index("-@", f"{THREADS}", str(path_bam))

    # Remove temporary files
    sam_temp = Path(TEMPDIR, "report", "bam").glob(f"tmp{randomnum}*.sam")
    [s.unlink() for s in sam_temp]
