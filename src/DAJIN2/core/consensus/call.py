from __future__ import annotations
from bisect import bisect_left
from collections import defaultdict
import midsv
from pathlib import Path
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


def transpose(cssplit):
    cs_transposed = []
    for cs in cssplit:
        cs_transposed.append(cs["CSSPLIT"].split(","))
    return list(zip(*cs_transposed))


def call_percentage(cssplit_sample: list[str], cssplit_control: list[str]) -> list[dict]:
    cssplit_control_transposed = transpose(cssplit_control)
    cssplit_sample_transposed = transpose(cssplit_sample)
    readnum_sample = len(cssplit_sample)
    readnum_control = len(cssplit_control)
    con_percentage = []
    for i in range(len(cssplit_sample_transposed)):
        d_control = defaultdict(int)
        for cs in cssplit_control_transposed[i]:
            d_control[cs] += 1 / readnum_control
        d_sample = defaultdict(int)
        for cs in cssplit_sample_transposed[i]:
            d_sample[cs] += 1 / readnum_sample
        d_sample_subtracted = defaultdict(int)
        for key, value in d_sample.items():
            if key.startswith("="):
                d_sample_subtracted[key] = value
                continue
            value = max(0, value - d_control[key])
            if value > 0:
                d_sample_subtracted[key] = value
        d_sum = sum(d_sample_subtracted.values())
        d = {key: value * 1 / d_sum for key, value in d_sample_subtracted.items()}
        d = dict(sorted(d.items(), key=lambda x: x[1], reverse=True))
        con_percentage.append(d)
    return con_percentage


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


def call_allele_name(keys: tuple, cons_seq: str, dict_allele: dict) -> str:
    ALLELE, SV, LABEL = keys
    allele_name = f"allele{LABEL}_{ALLELE}"
    if SV:
        allele_name += "_sv"
    elif cons_seq == dict_allele[ALLELE]:
        allele_name += "_intact"
    else:
        allele_name += "_mutated"
    return allele_name


def call(TEMPDIR, clust_sample, DICT_ALLELE, SAMPLE_NAME, CONTROL_NAME):
    cssplit_control = midsv.read_jsonl(Path(TEMPDIR, "midsv", f"{CONTROL_NAME}_control.jsonl"))
    cssplit_sample = midsv.read_jsonl(Path(TEMPDIR, "midsv", f"{SAMPLE_NAME}_control.jsonl"))
    cssplit_sample = join_listdicts(clust_sample, cssplit_sample, key="QNAME")
    cssplit_sample.sort(key=lambda x: (x["ALLELE"], x["SV"], x["LABEL"]))
    cons_percentage = defaultdict(list)
    cons_sequence = defaultdict(list)
    for keys, cssplits in groupby(cssplit_sample, key=lambda x: (x["ALLELE"], x["SV"], x["LABEL"])):
        cssplits = list(cssplits)
        cons_per = call_percentage(cssplits, cssplit_control)
        cons_seq = call_sequence(cons_per)
        allele_name = call_allele_name(keys, cons_seq, DICT_ALLELE)
        cons_percentage[allele_name] = cons_per
        cons_sequence[allele_name] = cons_seq
        for cs in cssplits:
            cs["NAME"] = allele_name
    for res in cssplit_sample:
        del res["RNAME"]
        del res["CSSPLIT"]
    cssplit_sample.sort(key=lambda x: x["LABEL"])
    return cssplit_sample, cons_percentage, cons_sequence

