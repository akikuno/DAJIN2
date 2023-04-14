from __future__ import annotations

import random
import re
from collections import Counter, defaultdict
from copy import deepcopy
from pathlib import Path
import midsv
from scipy import stats
from scipy.spatial.distance import cosine

"""
- 5-merでコントロールと比較して、コントロールにもあるプロファイルの場合はマッチで補正する
    - 断片のリードで、端から続くNは無視
    - 補正する前にすでにすべてがマッチなら即continue
"""


# def extract_indexes_with_both_ends_not_N(cssplits: list[list[str]]) -> list[tuple[int, int]]:
#     indexes = []
#     for cssplit in cssplits:
#         cssplit = ",".join(cssplit)
#         n_prefix = re.search(r"^(N,)+", cssplit)
#         left_idx = n_prefix.end() if n_prefix else 0
#         n_suffix = re.search(r"(,N)+$", cssplit)
#         right_idx = n_suffix.start() if n_suffix else len(cssplit)
#         # output index of splitted cssplit
#         left_idx = cssplit[:left_idx].count(",")
#         right_idx = cssplit.count(",") - cssplit[right_idx:].count(",")
#         indexes.append((left_idx, right_idx))
#     return indexes


# def call_count(cssplits: list[list[str]], indexes: list[tuple(int, int)]) -> dict[dict[str, int]]:
#     """Count cssplits within 3-mer range at mutation (or sequence error) loci.
#     """
#     count_kmer = defaultdict(Counter)
#     for cssplit, idx in zip(cssplits, indexes):
#         left_idx, right_idx = idx
#         for i in range(left_idx + 1, right_idx):
#             if cssplit[i].startswith("="):
#                 continue
#             kmer = ",".join([cssplit[i - 1], cssplit[i], cssplit[i + 1]])
#             count_kmer[i] += Counter([kmer])
#     counts = {i: dict(count_kmer[i]) for i in count_kmer.keys()}
#     return counts


# def call_percentage(cssplits: list[list[str]], counts: dict[dict[str, int]]) -> dict[dict[str, float]]:
#     coverage = len(cssplits)
#     percents = deepcopy(counts)
#     for i, c in counts.items():
#         for kmer, count in c.items():
#             percents[i][kmer] = count / coverage * 100
#     return percents


# def subtract_percentage(percent_sample: dict, percent_control: dict) -> dict[dict[str, float]]:
#     percent_subtracted = deepcopy(percent_sample)
#     for i, samp in percent_sample.items():
#         cont = percent_control.get(i)
#         for kmer_samp, per_samp in samp.items():
#             if cont.get(kmer_samp):
#                 percent_subtracted[i][kmer_samp] = per_samp - cont[kmer_samp]
#     return percent_subtracted


# def select_candidate_mutation(percent_subtracted: dict, threshold: float = 0.5) -> dict[dict[str, float]]:
#     candidate_mutation = dict()
#     for i, samp in percent_subtracted.items():
#         mutation = {k for k, v in samp.items() if v > threshold}
#         candidate_mutation.update({i: mutation})
#     return candidate_mutation


# def update_cssplits(cssplits: list, sequence: str, candidate_mutation: dict) -> list[list[str]]:
#     cssplits_update = deepcopy(cssplits)
#     for j, mutation in candidate_mutation.items():
#         for i, cssplit in enumerate(cssplits):
#             kmer = ",".join([cssplit[j - 1], cssplit[j], cssplit[j + 1]])
#             if kmer in mutation:
#                 continue
#             else:
#                 cssplits_update[i][j] = "=" + sequence[j]
#     return cssplits_update


###############################################################################
# Scratch... 5mers
###############################################################################


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


def count_indels_5mer(cssplits: list[list[str]], left_idx: int, right_idx: int) -> list[dict]:
    transposed = [list(t) for t in zip(*cssplits)]
    count_indels_5mer = []
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
        count_indels_5mer.append(count)
    return count_indels_5mer


def extractr_sequence_errors(count_5mer_sample, count_5mer_control):
    sequence_errors = [set() for _ in range(len(count_5mer_sample))]
    for i in range(len(sequence_errors)):
        for ids in ["ins", "del", "sub"]:
            x = count_5mer_sample[i][ids]
            y = count_5mer_control[i][ids]
            distance = 1 - cosine(x, y)
            _, pvalue = stats.ttest_ind(x, y, equal_var=False)
            if distance > 0.95 and pvalue > 0.05:
                sequence_errors[i].add(ids)
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
        for idx, cssplit in enumerate(cssplits):
            key = ",".join([cssplit[i - 1], cssplit[i + 1]])
            if cssplit[i] == "@":
                cssplits_atmark[idx] = key
            else:
                cssplits_sampling_key[key].append(cssplit[i])
                cssplits_sampling_all.append(cssplit[i])
        for idx, key in cssplits_atmark.items():
            if cssplits_sampling_key[key]:
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


def correct_sequence_error(TEMPDIR: Path, FASTA_ALLELES: dict[str, str], CONTROL_NAME: str, SAMPLE_NAME: str) -> None:
    for allele, sequence in FASTA_ALLELES.items():
        midsv_sample = midsv.read_jsonl((Path(TEMPDIR, "midsv", f"{SAMPLE_NAME}_splice_{allele}.jsonl")))
        midsv_control = midsv.read_jsonl((Path(TEMPDIR, "midsv", f"{CONTROL_NAME}_splice_{allele}.jsonl")))
        cssplits_control = [cs["CSSPLIT"].split(",") for cs in midsv_control]
        cssplits_sample = [cs["CSSPLIT"].split(",") for cs in midsv_sample]
        left_idx, right_idx = set_indexes(sequence)
        count_5mer_sample = count_indels_5mer(cssplits_sample, left_idx, right_idx)
        count_5mer_control = count_indels_5mer(cssplits_control, left_idx, right_idx)
        sequence_errors = extractr_sequence_errors(count_5mer_sample, count_5mer_control)
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


# allele = "control"
# sequence = FASTA_ALLELES[allele]
# midsv_sample = midsv.read_jsonl("DAJINResults/.tempdir/tyr-pm/midsv/barcode31_splice_control.jsonl")
# midsv_control = midsv.read_jsonl("DAJINResults/.tempdir/tyr-pm/midsv/barcode32_splice_control.jsonl")

# cssplits_sample = [cs["CSSPLIT"].split(",") for cs in midsv_sample]
# cssplits_control = [cs["CSSPLIT"].split(",") for cs in midsv_control]

# sequence = sequence[51 - 9 : 53]
# sequence_length = len(sequence)
# test_sample = []
# for cs in cssplits_sample:
#     test_sample.append(cs[51 - 9 : 53])

# test_control = []
# for cs in cssplits_control:
#     test_control.append(cs[51 - 9 : 53])


# count_5mer_sample = count_indels_5mer(test_sample, left_idx, right_idx)
# count_5mer_control = count_indels_5mer(test_control, left_idx, right_idx)

# sequence_errors = extractr_sequence_errors(count_5mer_sample, count_5mer_control)
# cssplits_sample_error_replaced = replace_errors_to_atmark(cssplits_sample, sequence_errors, left_idx, right_idx)
# cssplits_control_error_replaced = replace_errors_to_atmark(cssplits_control, sequence_errors, left_idx, right_idx)
# cssplits_error_replaced = replace_errors_to_atmark(test_sample, sequence_errors, left_idx, right_idx)


# cssplits_replaced = replace_at_loci_to_match(cssplits_replaced, sequence)


# cssplits = [["=T", "@", "=G"], ["=T", "@", "=G"], ["=T", "=T", "=G"], ["=T", "-T", "=G"]]

# cssplits_sample_atmark_replaced = replace_atmark(cssplits_sample_error_replaced, sequence)
# cssplits_control_atmark_replaced = replace_atmark(cssplits_control_error_replaced, sequence)

# before = defaultdict(int)
# for cs in cssplits_sample:
#     before[cs[51]] += 1

# before = sorted(before.items(), key=lambda x: -x[1])

# after = defaultdict(int)
# for cs in cssplits_atmark_replaced:
#     after[cs[51]] += 1

# after = sorted(after.items(), key=lambda x: -x[1])

# before
# after


# before = defaultdict(int)
# for cs in test_sample:
#     before[cs[-2]] += 1

# before = sorted(before.items(), key=lambda x: -x[1])

# after = defaultdict(int)
# for cs in cssplits_atmark_replaced:
#     after[cs[-2]] += 1

# after = sorted(after.items(), key=lambda x: -x[1])

# before
# after
# d = defaultdict(int)
# for cs in cssplits_sample:
#     d[cs[51]] += 1

# d
# cssplits_sample[0][51]

# before = defaultdict(int)
# for cs in cssplits_sample:
#     before[cs[828]] += 1

# before

# after = defaultdict(int)
# for cs in cssplits_sample_atmark_replaced:
#     after[cs[828]] += 1

# before
# after
