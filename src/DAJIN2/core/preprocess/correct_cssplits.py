from __future__ import annotations
from collections import Counter
import numpy as np
import re
from itertools import chain
from copy import deepcopy
from collections import defaultdict
from scipy.spatial.distance import cosine
from pathlib import Path

import midsv


def transpose(cssplits):
    return [list(cs) for cs in zip(*cssplits)]


def call_count(transpose_cssplits: list[list[str]]) -> list[dict[str:int]]:
    cssplit_counts = []
    for cssplit in transpose_cssplits:
        count = Counter(cssplit)
        count = dict(count)
        cssplit_counts.append(count)
    return cssplit_counts


def count_deletons_in_repaet(count_control, repeat_span):
    list_deletions = []
    count_deletions = []
    repeat_idx = 0
    repeat_start = repeat_span[repeat_idx][0]
    repeat_end = repeat_span[repeat_idx][1]
    for i, counts in enumerate(count_control):
        if repeat_start <= i < repeat_end:
            count_del = sum(val for key, val in counts.items() if key.startswith("-")) + 1
            count_deletions.append(count_del)
        elif i == repeat_end:
            list_deletions.append(count_deletions)
            repeat_idx += 1
            if repeat_idx == len(repeat_span):
                break
            repeat_start = repeat_span[repeat_idx][0]
            repeat_end = repeat_span[repeat_idx][1]
            count_deletions = []
            if repeat_start == i:
                count_del = sum(val for key, val in counts.items() if key.startswith("-"))
                count_deletions.append(count_del)
    return list_deletions


def find_repetitive_dels(count_control, count_sample, sequence) -> set[int]:
    repeat_regrex = r"A{4,}|C{4,}|G{4,}|T{4,}|N{4,}"
    repeat_span = list(x.span() for x in re.finditer(repeat_regrex, sequence))
    control_repdels = count_deletons_in_repaet(count_control, repeat_span)
    sample_repdels = count_deletons_in_repaet(count_sample, repeat_span)
    cossims = [1 - cosine(c, s) for c, s in zip(control_repdels, sample_repdels)]
    high_cossims = [list(range(span[0], span[1])) for cossim, span in zip(cossims, repeat_span) if cossim > 0.95]
    return set(chain.from_iterable(high_cossims))


def sampling(cnt: Counter, size: int):
    elements = []
    probs = []
    coverage = sum(cnt.values())
    for key, val in cnt.items():
        elements.append(key)
        probs.append(val / coverage)
    np.random.seed(1)
    samples = np.random.choice(a=elements, size=size, p=probs)
    return samples


def replace_repdels(transpose_cssplits: list[list[str]], repeat_dels: set):
    replased_cssplits = deepcopy(transpose_cssplits)
    for i, cssplits in enumerate(replased_cssplits):
        if i not in repeat_dels:
            continue
        cnt = Counter(cssplits)
        size = sum(1 for cs in cssplits if cs.startswith("-"))
        if size == 0:
            continue
        for key in list(cnt.keys()):
            if key.startswith("-"):
                del cnt[key]
        if cnt:
            samples = sampling(cnt, size)
        else:
            samples = ["N"] * size
        iter_samples = iter(samples)
        for j, cs in enumerate(cssplits):
            if cs.startswith("-"):
                replased_cssplits[i][j] = next(iter_samples)
    return replased_cssplits


# def replace_both_ends_n(transpose_cssplits: list[list[str]]):
#     d_samples = defaultdict(iter)
#     for i, cssplits in enumerate(transpose_cssplits):
#         cnt = Counter(cssplits)
#         if cnt["N"] == 0:  # No N
#             continue
#         if len(cnt) == 1 and cnt["N"] > 0:  # All N
#             continue
#         size = cnt["N"]
#         del cnt["N"]
#         samples = sampling(cnt, size)
#         d_samples[i] = iter(samples)
#     cssplits_replaced = [list(cs) for cs in zip(*transpose_cssplits)]
#     for i, cssplits in enumerate(cssplits_replaced):
#         for j, cs in enumerate(cssplits):
#             if cs != "N":
#                 break
#             cssplits_replaced[i][j] = next(d_samples[j])
#     for i, cssplits in enumerate(cssplits_replaced):
#         cssplits = cssplits[::-1]
#         for j, cs in enumerate(cssplits):
#             if cs != "N":
#                 break
#             try:
#                 cssplits_replaced[i][len(cssplits) - 1 - j] = next(d_samples[len(cssplits) - 1 - j])
#             except TypeError:
#                 pass
#     cssplits_replaced = [list(cs) for cs in zip(*cssplits_replaced)]
#     return cssplits_replaced


###############################################################################
# main
###############################################################################


def correct_cssplits(TEMPDIR: Path, FASTA_ALLELS: dict[str, str], CONTROL_NAME: str, SAMPLE_NAME: str) -> None:
    """
    - Mask the following mutations as sequence errors
        - 1. Repetitive deletions
        - 2. Unknown bases at the both end of sequence
        - 3. Common mutations between control and sample (<- TODO!)
    - Random sampling
    """
    for allele, sequence in FASTA_ALLELS.items():
        midsv_control = midsv.read_jsonl((Path(TEMPDIR, "midsv", f"{CONTROL_NAME}_{allele}.jsonl")))
        cssplits_control = [cs["CSSPLIT"].split(",") for cs in midsv_control]
        # Sample
        midsv_sample = midsv.read_jsonl((Path(TEMPDIR, "midsv", f"{SAMPLE_NAME}_{allele}.jsonl")))
        cssplits_sample = [cs["CSSPLIT"].split(",") for cs in midsv_sample]
        transpose_control = transpose(cssplits_control)
        transpose_sample = transpose(cssplits_sample)
        # Make count matrix
        count_control = call_count(transpose_control)
        count_sample = call_count(transpose_sample)
        # Find repetitive dels
        repeat_dels = find_repetitive_dels(count_control, count_sample, sequence)
        transpose_mask_control = replace_repdels(transpose_control, repeat_dels)
        transpose_mask_sample = replace_repdels(transpose_sample, repeat_dels)
        # transpose_mask_control = replace_both_ends_n(transpose_mask_control)
        # transpose_mask_sample = replace_both_ends_n(transpose_mask_sample)
        # Replace CSSPLIT
        transpose_mask_control = [",".join(t) for t in transpose(transpose_mask_control)]
        transpose_mask_sample = [",".join(t) for t in transpose(transpose_mask_sample)]
        for i, cssplits in enumerate(transpose_mask_control):
            midsv_control[i]["CSSPLIT"] = cssplits
        for i, cssplits in enumerate(transpose_mask_sample):
            midsv_sample[i]["CSSPLIT"] = cssplits
        # Save as a json
        if "QSCORE" in midsv_control[0]:
            for cont in midsv_control:
                del cont["QSCORE"]
        if "QSCORE" in midsv_sample[0]:
            for samp in midsv_sample:
                del samp["QSCORE"]
        midsv.write_jsonl(midsv_control, Path(TEMPDIR, "midsv", f"{CONTROL_NAME}_{allele}.jsonl"))
        midsv.write_jsonl(midsv_sample, Path(TEMPDIR, "midsv", f"{SAMPLE_NAME}_{allele}.jsonl"))