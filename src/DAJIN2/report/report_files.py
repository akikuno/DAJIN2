from __future__ import annotations
import textwrap
import cstag
from collections import defaultdict


def to_fasta(header: str, cons_seq: str) -> str:
    header = ">" + header
    cons_seq_wrap = textwrap.wrap(cons_seq, 80)
    fasta = "\n".join([header, *cons_seq_wrap]) + "\n"
    return fasta


def to_html(header: str, cons_per: list[dict]) -> str:
    cons_cssplit = [max(cons, key=cons.get) for cons in cons_per]
    prev_cstag = cons_cssplit[0]
    cons_cstag = []
    cons_cstag.append(prev_cstag)
    for current_cstag in cons_cssplit[1:]:
        if "=" == prev_cstag[0] == current_cstag[0]:
            cons_cstag.append(current_cstag[1])
        elif "-" == prev_cstag[0] == current_cstag[0]:
            cons_cstag.append(current_cstag[1])
        elif "+" == current_cstag[0]:
            cs_ins = current_cstag.replace("+", "").split("|")
            ins = "+" + "".join(cs_ins[:-1]) + cs_ins[-1]
            cons_cstag.append(ins)
        else:
            cons_cstag.append(current_cstag)
        prev_cstag = current_cstag
    cons_cstag = "".join(cons_cstag)
    return cstag.to_html(cons_cstag, header)


def to_vcf(header: str, cons_per: list[dict]) -> str:
    pass


def left_join(sam_contents: list[str], clust_sample: list[dict]) -> dict[list]:
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
            key = f'{read_right["ALLELE"]}-{read_right["SV"]}-{read_right["LABEL"]}'
            sam_groups[key].append(read_left)
            idx_left += 1
        else:
            idx_right += 1
    return sam_groups

