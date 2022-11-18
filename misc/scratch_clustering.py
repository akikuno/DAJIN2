from __future__ import annotations
from collections import Counter
import midsv
import numpy as np
import re
from pathlib import Path
from itertools import chain
from copy import deepcopy
import statsmodels.api as sm
import bisect
from collections import defaultdict
from pprint import pprint
from sklearn.cluster import MeanShift


np.random.seed(1)

# Stx2 deletion
# - deletion: 2012-2739 (727 bases)
# tests/data/knockout/test_barcode25.fq.gz
# tests/data/knockout/test_barcode30.fq.gz
# tests/data/knockout/design_stx2.fa


def call_count(transpose_cssplits: list[list[str]]) -> list[dict[str:int]]:
    cssplit_counts = []
    for cssplit in transpose_cssplits:
        count = Counter(cssplit)
        count = dict(count)
        cssplit_counts.append(count)
    return cssplit_counts


###############################################################################

allele = "deletion"
sv = True
sequence = DICT_ALLELE[allele]
knockin_loci = KNOCKIN_LOCI[allele]

# Control
midsv_control = midsv.read_jsonl((Path(TEMPDIR, "midsv", f"{CONTROL_NAME}_{allele}.jsonl")))
cssplits_control = [cs["CSSPLIT"].split(",") for cs in midsv_control]
transpose_control = list(zip(*cssplits_control))
count_control = call_count(transpose_control)

# Sample
cssplits_sample = [cs["CSSPLIT"].split(",") for cs in classif_sample if cs["ALLELE"] == allele and cs["SV"] == sv]
transpose_sample = list(zip(*cssplits_sample))
count_sample = call_count(transpose_sample)

# Find repetitive dels
repeat_dels = find_repetitive_dels(count_control, count_sample, sequence)  # *依存


qscores = [cs["QSCORE"].split(",") for cs in midsv_control]

cssplits_replaced = replace_n_to_match(cssplits_control, sequence)
cssplits_replaced = replace_lowquality_to_match(cssplits_replaced, qscores, sequence)
transpose_control = transpose(cssplits_replaced)
count_control = call_count(transpose_control)

knockin_loci = find_knockin_loci(count_control, sequence)
count_control = recount_knockin_loci(count_control, knockin_loci, sequence)

# Sample
cssplits_sample = [cs["CSSPLIT"].split(",") for cs in classif_sample if cs["ALLELE"] == allele and cs["SV"] == sv]
transpose_sample = transpose(cssplits_sample)
qscores = [cs["QSCORE"].split(",") for cs in classif_sample if cs["ALLELE"] == allele and cs["SV"] == sv]
cssplits_replaced = replace_n_to_match(cssplits_sample, sequence)
cssplits_replaced = replace_lowquality_to_match(cssplits_replaced, qscores, sequence)
transpose_sample = transpose(cssplits_replaced)
count_sample = call_count(transpose_sample)
