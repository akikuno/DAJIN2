from __future__ import annotations

# from bisect import bisect_left
from collections import defaultdict
from itertools import groupby

# def merge_dict1(d1, d2):
#     d_merged = d1.copy()
#     for k, v in d2.items():
#         if k not in d_merged:
#             d_merged[k] = v
#     return d_merged


# def join_listdicts(list_dict1: list[dict], list_dict2: list[dict], key: str) -> list[dict]:
#     list_dict2.sort(key=lambda x: x[key])
#     key_dict2 = [d[key] for d in list_dict2]
#     dict_merged = []
#     for d1 in list_dict1:
#         key_dict1 = d1[key]
#         if key_dict1 not in key_dict2:
#             continue
#         idx = bisect_left(key_dict2, key_dict1)
#         d2 = list_dict2[idx]
#         dict_merged.append(merge_dict1(d1, d2))
#     return dict_merged


def transpose(cssplits: list[list[str]]) -> list[list[str]]:
    return [list(cs) for cs in zip(*cssplits)]


def call_percentage(cssplits: list[str]) -> list[dict[str, float]]:
    """
    Call position weight matrix in defferent loci.
    Non defferent loci are annotated to "Match" or "Unknown(N)"
    """
    cssplits_transposed = transpose(cssplits)
    coverage = len(cssplits)
    cons_percentage = []
    for cs_transposed in cssplits_transposed:
        count_cs = defaultdict(int)
        for cs in cs_transposed:
            count_cs[cs] += 1 / coverage * 100
        count_cs_sorted = dict(sorted(count_cs.items(), key=lambda x: x[1], reverse=True))
        cons_percentage.append(count_cs_sorted)
    return cons_percentage


def call_sequence(cons_percentage: list[dict[str, float]]) -> str:
    consensus_sequence = []
    for cons_per in cons_percentage:
        cons = max(cons_per, key=cons_per.get)
        if cons.startswith("="):
            cons = cons.replace("=", "")
        elif cons.startswith("-"):
            continue
        elif cons.startswith("*"):
            cons = cons[-1]
        elif cons.startswith("+"):
            cons_ins = cons.split("|")
            if cons_ins[-1].startswith("="):
                cons = cons.replace("=", "")
            elif cons_ins[-1].startswith("-"):
                cons = "".join(cons_ins[:-1])
            elif cons_ins[-1].startswith("*"):
                cons = "".join([*cons_ins[:-1], cons_ins[-1][-1]])
            cons = cons.replace("+", "")
            cons = cons.replace("|", "")
        consensus_sequence.append(cons)
    return "".join(consensus_sequence)


def call_allele_name(keys: tuple, cons_seq: str, FASTA_ALLELES: dict) -> str:
    ALLELE, SV, LABEL, PERCENT = keys
    allele_name = f"allele{LABEL}_{ALLELE}"
    if SV:
        allele_name += "_sv"
    elif cons_seq == FASTA_ALLELES[ALLELE]:
        allele_name += "_intact"
    else:
        allele_name += "_mutated"
    allele_name += f"_{PERCENT}%"
    return allele_name


def call(clust_sample: list[dict], FASTA_ALLELES: dict):
    """
    - RESULT_SAMPLE: a list of dicts including "QNAME", "ALLELE", "CSSPLIT", "SV", "LABEL", "READNUM", "PERCENT", and "NAME"
    """
    RESULT_SAMPLE = []
    cons_percentage = defaultdict(list)
    cons_sequence = defaultdict(list)
    clust_sample.sort(key=lambda x: x["LABEL"])
    for _, group in groupby(clust_sample, key=lambda x: x["LABEL"]):
        clust = list(group)
        keys = (
            clust[0]["ALLELE"],
            clust[0]["SV"],
            clust[0]["LABEL"],
            clust[0]["PERCENT"],
        )
        cssplits = [cs["CSSPLIT"].split(",") for cs in clust]
        cons_per = call_percentage(cssplits)
        cons_seq = call_sequence(cons_per)
        allele_name = call_allele_name(keys, cons_seq, FASTA_ALLELES)
        cons_percentage[allele_name] = cons_per
        cons_sequence[allele_name] = cons_seq
        for cs in clust:
            cs["NAME"] = allele_name
        RESULT_SAMPLE.extend(clust)
    RESULT_SAMPLE.sort(key=lambda x: x["LABEL"])
    return RESULT_SAMPLE, cons_percentage, cons_sequence
