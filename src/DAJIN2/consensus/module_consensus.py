from __future__ import annotations
import textwrap
from collections import Counter
from bisect import bisect_left
from collections import defaultdict


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


def percentage(cssplit_sample: list[str], cssplit_control: list[str]) -> list[dict]:
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


def call(cssplits: list[str], diffloci: list[int]) -> list[str]:
    cs = [cs.split(",") for cs in cssplits]
    cons = []
    for difflocus in diffloci:
        cons_difflocus = []
        for c in cs:
            cons_difflocus.append(c[difflocus])
        cons.append(Counter(cons_difflocus).most_common()[0][0])
    return cons


def call_fasta_seq(cons_percentage_by_key: list[dict]) -> str:
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


def call_fasta(key: str, cons_percentage_by_key: list[dict]) -> str:
    header = ">" + key
    cons_seq = call_fasta_seq(cons_percentage_by_key)
    cons_seq_wrap = textwrap.wrap(cons_seq, 80)
    fasta = "\n".join([header, *cons_seq_wrap]) + "\n"
    return fasta
