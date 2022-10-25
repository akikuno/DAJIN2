from __future__ import annotations
from bisect import bisect_left
from collections import defaultdict
from itertools import groupby


def merge_dict1(d1, d2):
    d_merged = d1.copy()
    for k, v in d2.items():
        if k not in d_merged:
            d_merged[k] = v
    return d_merged


def join_listdicts(list_dict1: list[dict], list_dict2: list[dict], key: str) -> list[dict]:
    list_dict2.sort(key=lambda x: x[key])
    key_dict2 = [d[key] for d in list_dict2]
    dict_merged = []
    for d1 in list_dict1:
        key_dict1 = d1[key]
        if key_dict1 not in key_dict2:
            continue
        idx = bisect_left(key_dict2, key_dict1)
        d2 = list_dict2[idx]
        dict_merged.append(merge_dict1(d1, d2))
    return dict_merged


def transpose(cssplit: list[str]) -> list[tuple[str]]:
    cs_transposed = []
    for cs in cssplit:
        cs_transposed.append(cs["CSSPLIT"].split(","))
    return list(zip(*cs_transposed))


def call_percentage(cssplit_sample: list[str], diffloci: list[int], repetitive_del_loci: list[int]) -> list[dict]:
    """
    Call position weight matrix in defferent loci.
    Non defferent loci are annotated to "Match" or "Unknown(N)"
    """
    cssplit_sample_transposed = transpose(cssplit_sample)
    coverage = len(cssplit_sample)
    cons_percentage = []
    for idx, cs_transposed in enumerate(cssplit_sample_transposed):
        if idx in diffloci:
            count_cs = defaultdict(int)
            for cs in cs_transposed:
                if cs.startswith("-") and idx in repetitive_del_loci:
                    continue
                count_cs[cs] += 1 / coverage
            if not count_cs:
                count_cs["N"] = 1.0
            count_cs_sorted = dict(sorted(count_cs.items(), key=lambda x: x[1], reverse=True))
            cons_percentage.append(count_cs_sorted)
        else:
            # return Match when one match found.
            count_cs = {"N": 1.0}
            for cs in cs_transposed:
                if cs.startswith("="):
                    count_cs = {cs: 1.0}
                    break
            cons_percentage.append(count_cs)
    return cons_percentage


def call_sequence(cons_percentage_by_key: list[dict]) -> str:
    consensus_sequence = []
    for cons_per in cons_percentage_by_key:
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


def call_allele_name(keys: tuple, cons_seq: str, DICT_ALLELE: dict) -> str:
    ALLELE, SV, LABEL, PERCENT = keys
    allele_name = f"allele{LABEL}_{ALLELE}"
    if SV:
        allele_name += "_sv"
    elif cons_seq == DICT_ALLELE[ALLELE]:
        allele_name += "_intact"
    else:
        allele_name += "_mutated"
    allele_name += f"_{PERCENT}%"
    return allele_name


def call(
    clust_sample: list[dict], DIFFLOCI_ALLELES: defaultdict[int], REPETITIVE_DELLOCI: list(tuple), DICT_ALLELE: dict
):
    result_sample = []
    cons_percentage = defaultdict(list)
    cons_sequence = defaultdict(list)
    clust_sample.sort(key=lambda x: x["LABEL"])
    for _, group in groupby(clust_sample, key=lambda x: x["LABEL"]):
        cssplit_sample = list(group)
        keys = (
            cssplit_sample[0]["ALLELE"],
            cssplit_sample[0]["SV"],
            cssplit_sample[0]["LABEL"],
            cssplit_sample[0]["PERCENT"],
        )
        diffloci = DIFFLOCI_ALLELES[f"{keys[0]}-{keys[1]}"]
        repetitive_del_loci = REPETITIVE_DELLOCI[f"{keys[0]}-{keys[1]}"]
        cons_per = call_percentage(cssplit_sample, diffloci, repetitive_del_loci)
        cons_seq = call_sequence(cons_per)
        allele_name = call_allele_name(keys, cons_seq, DICT_ALLELE)
        cons_percentage[allele_name] = cons_per
        cons_sequence[allele_name] = cons_seq
        for cs in cssplit_sample:
            cs["NAME"] = allele_name
        result_sample.extend(cssplit_sample)
    result_sample.sort(key=lambda x: x["LABEL"])
    return result_sample, cons_percentage, cons_sequence


# def call(TEMPDIR, clust_sample, DICT_ALLELE, SAMPLE_NAME, CONTROL_NAME):
#     result_sample = []
#     cons_percentage = defaultdict(list)
#     cons_sequence = defaultdict(list)
#     clust_sample.sort(key=lambda x: x["LABEL"])
#     for _, group in groupby(clust_sample, key=lambda x: x["LABEL"]):
#         clust = list(group)
#         keys = (clust[0]["ALLELE"], clust[0]["SV"], clust[0]["LABEL"], clust[0]["PERCENT"])
#         allele = keys[0]
#         cssplit_control = midsv.read_jsonl(Path(TEMPDIR, "midsv", f"{CONTROL_NAME}_{allele}.jsonl"))
#         cssplit_sample = midsv.read_jsonl(Path(TEMPDIR, "midsv", f"{SAMPLE_NAME}_{allele}.jsonl"))
#         cssplit_sample = join_listdicts(clust, cssplit_sample, key="QNAME")
#         cons_per = call_percentage(cssplit_sample, cssplit_control)
#         cons_seq = call_sequence(cons_per)
#         allele_name = call_allele_name(keys, cons_seq, DICT_ALLELE)
#         cons_percentage[allele_name] = cons_per
#         cons_sequence[allele_name] = cons_seq
#         for cs in cssplit_sample:
#             cs["NAME"] = allele_name
#             del cs["RNAME"]
#             del cs["CSSPLIT"]
#         result_sample.extend(cssplit_sample)
#     result_sample.sort(key=lambda x: x["LABEL"])
#     return result_sample, cons_percentage, cons_sequence

