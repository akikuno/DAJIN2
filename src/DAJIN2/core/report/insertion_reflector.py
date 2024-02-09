from __future__ import annotations

import cstag


def get_index_of_insertions(ref_cstag: str) -> set:
    index_of_insertions = set()
    start = 0
    for cs in cstag.split(ref_cstag):
        if cs.startswith("+"):
            index_of_insertions |= set(range(start, start + len(cs[1:])))
            start += len(cs[1:])
        elif cs.startswith("="):
            start += len(cs[1:])
        elif cs.startswith("*"):
            start += 1
        elif cs.startswith("-"):
            continue
    return index_of_insertions


def split_cstag(cs_tag: str) -> list[str]:
    cstag_split = []
    for cs in cstag.split(cs_tag):
        if cs.startswith("=") or cs.startswith("-") or cs.startswith("+"):
            op = cs[0]
            cstag_split += [op + c for c in list(cs[1:])]
        elif cs.startswith("*"):
            cstag_split.append(cs)
    return cstag_split


def apply_insertion(cs_split: list[str], index_of_insertions: set) -> list[str]:
    idx = -1
    cs_insertion = []
    for cs in cs_split:
        if not cs.startswith("+"):
            idx += 1
        if idx in index_of_insertions:
            if cs.startswith("-"):
                continue
            cs_insertion.append("+" + cs[-1].lower())
        else:
            cs_insertion.append(cs)
    return cs_insertion


def convert_to_cstag(cs_insertion: list[str]) -> str:
    cs_tag = []
    prev_op = "@"
    for cs in cs_insertion:
        if cs.startswith(prev_op) and prev_op != "*":
            cs_tag.append(cs[1:])
            continue
        prev_op = cs[0]
        cs_tag.append(cs)
    return "".join(cs_tag)


def reflect_ref_insertion_to_query(ref_cstag: str, que_cstag: str) -> str:
    """
    Translate the cstag of the query corresponding to the insertion site of the reference into an insertion.

    Example:
        ref_cstag = "=AC+aaa=GT"
        que_cstag = "=ACAAAGT"
        reflect_ref_insertion_to_query(ref_cstag, que_cstag)  # -> "=AC+aaa=GT"
    """
    index_of_insertions = get_index_of_insertions(ref_cstag)
    cs_split = split_cstag(que_cstag)
    cs_insertion = apply_insertion(cs_split, index_of_insertions)
    return convert_to_cstag(cs_insertion)
