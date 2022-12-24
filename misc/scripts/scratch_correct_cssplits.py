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

"""
- 5-merでコントロールと比較して、コントロールにもあるプロファイルの場合はマッチで補正する
    - 補正する前にすでにすべてがマッチなら即continue
    - 断片のリードで、端から続くNは無視
"""

