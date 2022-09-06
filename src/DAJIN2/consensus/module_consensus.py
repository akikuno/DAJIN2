from __future__ import annotations
from collections import Counter
from bisect import bisect_left


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


def call(cssplits: list[str], diffloci: list[int]) -> list[str]:
    cs = [cs.split(",") for cs in cssplits]
    cons = []
    for difflocus in diffloci:
        cons_difflocus = []
        for c in cs:
            cons_difflocus.append(c[difflocus])
        cons.append(Counter(cons_difflocus).most_common()[0][0])
    return cons
