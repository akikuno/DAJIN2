from __future__ import annotations
from collections import Counter
import midsv
import numpy as np
import re
from itertools import chain
from copy import deepcopy


#!--------------------------------------------------
print(f"mask_repeditive_dels.py: WOURKING IN PROGRESS...")
#!--------------------------------------------------


def pearson_corr(x: list[int], y: list[int]):
    x_diff = np.array(x) - np.mean(x)
    y_diff = np.array(y) - np.mean(y)
    if (np.sqrt(sum(x_diff ** 2)) * np.sqrt(sum(y_diff ** 2))) == 0:
        return 0
    return np.dot(x_diff, y_diff) / (np.sqrt(sum(x_diff ** 2)) * np.sqrt(sum(y_diff ** 2)))


def replaceNtoMatch(cssplit: list[str], sequence) -> list[str]:
    cssplit_replaced = deepcopy(cssplit)
    # Replace N to @ at the left ends
    for i, cs in enumerate(cssplit_replaced):
        if cs != "N":
            break
        cssplit_replaced[i] = "=" + sequence[i]
    # Replace N to @ at the right ends
    cssplit_replaced = cssplit_replaced[::-1]
    for i, cs in enumerate(cssplit_replaced):
        if cs != "N":
            break
        cssplit_replaced[i] = "=" + sequence[::-1][i]
    cssplit_replaced = cssplit_replaced[::-1]
    return cssplit_replaced


def call_count(cssplit_transposed):
    cssplit_counts = []
    for cssplit in cssplit_transposed:
        count = Counter(cssplit)
        count = dict(count)
        cssplit_counts.append(count)
    return cssplit_counts


def count_deletons_in_repaet(control_transposed, sequence):
    repeat_regrex = r"A{4,}|C{4,}|G{4,}|T{4,}|N{4,}"
    repeat_span = list(x.span() for x in re.finditer(repeat_regrex, sequence))
    repeats = [list(range(span[0], span[1])) for span in repeat_span]
    repeats_delcounts = []
    for reps in repeats:
        repeat_dels = []
        for rep in reps:
            repeat_dels.append(sum([1 for c in control_transposed[rep] if c.startswith("-")]))
        repeats_delcounts.append(repeat_dels)
    return repeats_delcounts


def find_repetitive_dels(control_count, sample_count, sequence):
    control_repdels = count_deletons_in_repaet(control_count, sequence)
    sample_repdels = count_deletons_in_repaet(sample_count, sequence)
    repeat_regrex = r"A{4,}|C{4,}|G{4,}|T{4,}|N{4,}"
    repeat_span = list(x.span() for x in re.finditer(repeat_regrex, sequence))
    corrs = [pearson_corr(c, s) for c, s in zip(control_repdels, sample_repdels)]
    x = [list(range(span[0], span[1])) for corr, span in zip(corrs, repeat_span) if corr > 0.95]
    return set(chain.from_iterable(x))


sequence = DICT_ALLELE["control"]
control_cssplit = midsv.read_jsonl(Path(TEMPDIR, "midsv", f"{CONTROL_NAME}_control.jsonl"))
coverage = len(control_cssplit)
control_cssplit = [cs["CSSPLIT"].split(",") for cs in control_cssplit]
control_cssplit = [replaceNtoMatch(cs, sequence) for cs in control_cssplit]
control_transposed = list(zip(*control_cssplit))
control_repdels = count_deletons_in_repaet(control_transposed, sequence)
average_left = sum([c[0] for c in control_repdels]) / len(control_repdels)
average_right = sum([c[-1] for c in control_repdels]) / len(control_repdels)
if average_left > average_right:
    control_repdels = [c[::-1] for c in control_repdels]

control_repdels_percent = []
for cont in control_repdels:
    control_repdels_percent.append([c / coverage for c in cont])

from itertools import zip_longest

control_repdels_transposed = list(zip_longest(*control_repdels))
control_repdels_statistics = []
for cont in control_repdels_transposed:
    cont_percent = [c / coverage * 100 for c in cont if c]
    cont_percent_quatile = np.quantile(cont_percent, 0.9)
    cont_percent_95 = np.std(cont_percent) * 2
    cont_percent_quatile += cont_percent_95
    control_repdels_statistics.append(cont_percent_quatile)

