from __future__ import annotations
from collections import Counter
import numpy as np
import re
from itertools import chain
from copy import deepcopy
from collections import defaultdict

from scipy import stats
from scipy.spatial.distance import cosine

from pathlib import Path

import midsv

"""
- 5-merでコントロールと比較して、コントロールにもあるプロファイルの場合はマッチで補正する
    - 断片のリードで、端から続くNは無視
    - 補正する前にすでにすべてがマッチなら即continue
"""


def extract_indexes_with_both_ends_not_N(cssplits: str) -> tuple(int, int):
    n_prefix = re.search(r"^(N,)+", cssplits)
    left_idx = n_prefix.end() if n_prefix else 0
    n_suffix = re.search(r"(,N)+$", cssplits)
    right_idx = n_suffix.start() if n_suffix else len(cssplits)
    # output index of splitted cssplits
    left_idx = cssplits[:left_idx].count(",")
    right_idx = cssplits.count(",") - cssplits[right_idx:].count(",")
    return left_idx, right_idx


# cssplits = "N,N,N,N,=A,=C,=G,=T"
# cssplits = "N,N,=A,=C,=G,=T,N,N"
# left_idx, right_idx = extract_indexes_with_both_ends_not_N(cssplits)
# # cssplits.split(",")[left_idx]
# # cssplits.split(",")[right_idx]

# cssplits_sample = ["=A,=C,=G,=T", "*AT,-C,+A|G,=t"]

# a = [0.01, 0.02, 0.03]
# b = [0.97, 0.98, 0.99]
# 1 - cosine(a, b)

# _, pval = stats.ttest_ind(a, b, equal_var=False)
# pval
