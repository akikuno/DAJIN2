from __future__ import annotations

import random
import re
from collections import Counter, defaultdict
from copy import deepcopy
from pathlib import Path
import midsv
import numpy as np
from scipy import stats
from scipy.spatial import distance


def set_indexes(sequence: str):
    sequence_length = len(sequence)
    num_subset = sequence_length % 5
    left_idx = 0
    right_idx = sequence_length
    if num_subset == 1:
        left_idx += 1
    elif num_subset == 2:
        left_idx += 1
        right_idx -= 1
    elif num_subset == 3:
        left_idx += 2
        right_idx -= 1
    elif num_subset == 4:
        left_idx += 2
        right_idx -= 2
    return left_idx, right_idx


def count_5mer_indels(cssplits: list[list[str]], left_idx: int, right_idx: int) -> list[dict]:
    transposed = [list(t) for t in zip(*cssplits)]
    count_5mer = []
    for i in range(left_idx, right_idx, 5):
        count = {"ins": [1] * 5, "del": [1] * 5, "sub": [1] * 5}
        cssplits_5mer = transposed[i : i + 5]
        for j, cs in enumerate(cssplits_5mer):
            counter = Counter(cs)
            for key, cnt in counter.items():
                if key.startswith("=") or key == "N" or re.search(r"a|c|g|t|n", key):
                    continue
                if key.startswith("+"):
                    count["ins"][j] += cnt
                elif key.startswith("-"):
                    count["del"][j] += cnt
                elif key.startswith("*"):
                    count["sub"][j] += cnt
        count_5mer.append(count)
    return count_5mer


def remove_minor_indels(cssplits: list[list[str]], count_5mer: list[dict]) -> list[dict]:
    coverage = len(cssplits)
    count_5mer_filtered = []
    for count in count_5mer:
        dict_mutation = defaultdict(list)
        for mutation in ["ins", "del", "sub"]:
            if all(True if c < coverage * 0.01 else False for c in count[mutation]):
                dict_mutation[mutation] = [1] * 5
            else:
                dict_mutation[mutation] = count[mutation]
        count_5mer_filtered.append(dict_mutation)
    return count_5mer_filtered


def extract_sequence_errors(count_5mer_sample, count_5mer_control):
    sequence_errors = [set() for _ in range(len(count_5mer_sample))]
    dists = defaultdict(list)
    # Calculate Jensen-Shannon distance
    for samp, cont in zip(count_5mer_sample, count_5mer_control):
        for mutation in ["ins", "del", "sub"]:
            s = samp[mutation]
            c = cont[mutation]
            dists[mutation].append(distance.jensenshannon(s, c))
    # Discrimitate seq errors and real mutation using Hotelling's T-squared distribution
    dists_all = np.array(list(dists.values())).flatten()
    avg = np.average(dists_all[~np.isnan(dists_all)])
    var = np.var(dists_all[~np.isnan(dists_all)])
    threshold = 0.05
    for mutation in ["ins", "del", "sub"]:
        dists_subset = dists[mutation]
        scores = [(xi - avg) ** 2 / var for xi in dists_subset]
        thres = stats.chi2.interval(1 - threshold, 1)[1]
        for i, score in enumerate(scores):
            # 'nan' means the two distributions have too different, so it could be a real mutation
            if np.isnan(score):
                continue
            if score < thres:
                sequence_errors[i].add(mutation)
    return sequence_errors


def replace_errors_to_atmark(cssplits_sample, sequence_errors, left_idx, right_idx):
    cssplits_replaced = []
    for samp in cssplits_sample:
        samp_replaced = deepcopy(samp)
        for idx_error, idx_5mer in enumerate(range(left_idx, right_idx, 5)):
            samp_5mer = samp[idx_5mer : idx_5mer + 5]
            error = sequence_errors[idx_error]
            if "ins" in error:
                samp_5mer = ["@" if cs.startswith("+") else cs for cs in samp_5mer]
            if "del" in error:
                samp_5mer = ["@" if cs.startswith("-") else cs for cs in samp_5mer]
            if "sub" in error:
                samp_5mer = ["@" if cs.startswith("*") else cs for cs in samp_5mer]
            samp_replaced[idx_5mer : idx_5mer + 5] = samp_5mer
        cssplits_replaced.append(samp_replaced)
    return cssplits_replaced


def replace_atmark(cssplits: list[list[str]], sequence: str) -> list[list[str]]:
    random.seed(1)
    cssplits_replaced = deepcopy(cssplits)
    sequence_length = len(sequence)
    for i in range(1, sequence_length - 1):
        cssplits_atmark = defaultdict(str)
        cssplits_sampling_key = defaultdict(list)
        cssplits_sampling_all = []
        flag_all_atmark = True
        for idx, cssplit in enumerate(cssplits):
            key = ",".join([cssplit[i - 1], cssplit[i + 1]])
            if cssplit[i] == "@":
                cssplits_atmark[idx] = key
            else:
                cssplits_sampling_key[key].append(cssplit[i])
                cssplits_sampling_all.append(cssplit[i])
                flag_all_atmark = False
        for idx, key in cssplits_atmark.items():
            if flag_all_atmark:
                cssplits_replaced[idx][i] = "N"
            elif cssplits_sampling_key[key]:
                cssplits_replaced[idx][i] = random.choice(cssplits_sampling_key[key])
            else:
                cssplits_replaced[idx][i] = random.choice(cssplits_sampling_all)
    for cs in cssplits_replaced:
        if cs[0] == "@":
            cs[0] = "=" + sequence[0]
        if cs[-1] == "@":
            cs[-1] = "=" + sequence[-1]
    return cssplits_replaced


###############################################################################
# main
###############################################################################


def execute(TEMPDIR: Path, FASTA_ALLELES: dict[str, str], CONTROL_NAME: str, SAMPLE_NAME: str) -> None:
    for allele, sequence in FASTA_ALLELES.items():
        midsv_sample = midsv.read_jsonl((Path(TEMPDIR, "midsv", f"{SAMPLE_NAME}_splice_{allele}.jsonl")))
        midsv_control = midsv.read_jsonl((Path(TEMPDIR, "midsv", f"{CONTROL_NAME}_splice_{allele}.jsonl")))
        cssplits_sample = [cs["CSSPLIT"].split(",") for cs in midsv_sample]
        cssplits_control = [cs["CSSPLIT"].split(",") for cs in midsv_control]
        # Extract sequence errors
        left_idx, right_idx = set_indexes(sequence)
        count_5mer_sample = count_5mer_indels(cssplits_sample, left_idx, right_idx)
        count_5mer_control = count_5mer_indels(cssplits_control, left_idx, right_idx)
        count_5mer_sample = remove_minor_indels(cssplits_sample, count_5mer_sample)
        count_5mer_control = remove_minor_indels(cssplits_control, count_5mer_control)
        sequence_errors = extract_sequence_errors(count_5mer_sample, count_5mer_control)
        # Correct sequence errors
        cssplits_sample_error_replaced = replace_errors_to_atmark(cssplits_sample, sequence_errors, left_idx, right_idx)
        cssplits_control_error_replaced = replace_errors_to_atmark(
            cssplits_control, sequence_errors, left_idx, right_idx
        )
        cssplits_sample_atmark_replaced = replace_atmark(cssplits_sample_error_replaced, sequence)
        cssplits_control_atmark_replaced = replace_atmark(cssplits_control_error_replaced, sequence)
        # Replace CSSPLIT
        cssplits_sample_corrected = [",".join(cs) for cs in cssplits_sample_atmark_replaced]
        cssplits_control_corrected = [",".join(cs) for cs in cssplits_control_atmark_replaced]
        for i, cssplits in enumerate(cssplits_sample_corrected):
            midsv_sample[i]["CSSPLIT"] = cssplits
        for i, cssplits in enumerate(cssplits_control_corrected):
            midsv_control[i]["CSSPLIT"] = cssplits
        midsv.write_jsonl(midsv_control, Path(TEMPDIR, "midsv", f"{CONTROL_NAME}_splice_{allele}.jsonl"))
        midsv.write_jsonl(midsv_sample, Path(TEMPDIR, "midsv", f"{SAMPLE_NAME}_splice_{allele}.jsonl"))
