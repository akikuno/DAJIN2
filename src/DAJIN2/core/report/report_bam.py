from __future__ import annotations

import random
from collections import defaultdict
from itertools import groupby
from pathlib import Path

import midsv
import pysam

from DAJIN2.utils import sam_handler


def realign(sam: list[list[str]], GENOME_COODINATES: dict) -> list[str]:
    sam_headers = [s for s in sam if s[0].startswith("@")]
    sam_contents = [s for s in sam if not s[0].startswith("@")]
    for s in sam_headers:
        if s[0] != "@SQ":
            continue
        s[1] = f'SN:{GENOME_COODINATES["chrom"]}'
        s[2] = f'LN:{GENOME_COODINATES["chrom_size"]}'
    for s in sam_contents:
        s[2] = GENOME_COODINATES["chrom"]
    if GENOME_COODINATES["strand"] == "-":
        sam_contents = sam_handler.revcomp_sam(sam_contents, GENOME_COODINATES["end"])
    else:
        for s in sam_contents:
            s[3] = str(int(s[3]) + GENOME_COODINATES["start"] - 1)
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


def write_sam_to_bam(sam: list[list[str]], path_sam: str | Path, path_bam: str | Path, threads: int = 1) -> None:
    formatted_sam = "\n".join("\t".join(s) for s in sam)
    Path(path_sam).write_text(formatted_sam + "\n")
    pysam.sort("-@", f"{threads}", "-o", str(path_bam), str(path_sam))
    pysam.index("-@", f"{threads}", str(path_bam))


def update_sam(sam: list, GENOME_COODINATES: dict = None) -> list:
    sam_update = sam.copy()
    sam_update = sam_handler.remove_overlapped_reads(sam_update)
    sam_update = sam_handler.remove_microhomology(sam_update)
    if GENOME_COODINATES["genome"]:
        sam_update = realign(sam_update, GENOME_COODINATES)
    return sam_update


def output_bam(TEMPDIR, NAME, GENOME_COODINATES, THREADS, RESULT_SAMPLE=None, is_control=False) -> None:
    randomnum = random.randint(100_000, 999_999)
    path_sam_input = Path(TEMPDIR, NAME, "sam", "map-ont_control.sam")
    sam = list(midsv.read_sam(path_sam_input))
    # Update sam
    sam_update = update_sam(sam, GENOME_COODINATES)
    # Output SAM and BAM
    path_sam_output = Path(TEMPDIR, "report", "BAM", f"tmp{randomnum}_{NAME}_control.sam")
    path_bam_output = Path(TEMPDIR, "report", "BAM", NAME, f"{NAME}.bam")
    write_sam_to_bam(sam_update, path_sam_output, path_bam_output, THREADS)
    # Prepare SAM headers and contents
    sam_headers = [s for s in sam_update if s[0].startswith("@")]
    sam_contents = [s for s in sam_update if not s[0].startswith("@")]
    if is_control:
        qnames = set(list(set(s[0] for s in sam_contents[:10000]))[:100])
        sam_subset = [s for s in sam_update if s[0] in qnames]
        path_sam_output = Path(TEMPDIR, "report", "BAM", f"tmp{randomnum}_{NAME}_control_cache.sam")
        path_bam_output = Path(TEMPDIR, "cache", ".igvjs", NAME, "control.bam")
        write_sam_to_bam(sam_headers + sam_subset, path_sam_output, path_bam_output, THREADS)
    else:
        sam_groups = group_by_name(sam_contents, RESULT_SAMPLE)
        qnames_by_name = subset_qnames(RESULT_SAMPLE)
        # Output SAM and BAM
        for name, sam_content in sam_groups.items():
            # BAM
            path_sam_output = Path(TEMPDIR, "report", "bam", f"tmp{randomnum}_{name}.sam")
            path_bam_output = Path(TEMPDIR, "report", "BAM", NAME, f"{NAME}_{name}.bam")
            write_sam_to_bam(sam_headers + sam_content, path_sam_output, path_bam_output, THREADS)
            # igvjs
            sam_subset = subset_reads(name, sam_content, qnames_by_name)
            path_sam_output = Path(TEMPDIR, "report", "bam", f"tmp{randomnum}_{name}_subset.sam")
            path_bam_output = Path(TEMPDIR, "report", ".igvjs", NAME, f"{name}.bam")
            write_sam_to_bam(sam_headers + sam_subset, path_sam_output, path_bam_output, THREADS)
    # Remove temporary files
    sam_temp = Path(TEMPDIR, "report", "BAM").glob(f"tmp{randomnum}*.sam")
    [s.unlink() for s in sam_temp]
