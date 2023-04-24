from __future__ import annotations

import random
import re
from collections import Counter
from collections import defaultdict
from pathlib import Path

import midsv
import numpy as np
from scipy import stats
from scipy.spatial import distance
from sklearn.mixture import GaussianMixture
from sklearn.neighbors import LocalOutlierFactor


def _count_indels(cssplits: list[list[str]]) -> dict[str, list[int]]:
    count = {"ins": [1] * len(cssplits[0]),
            "del": [1] * len(cssplits[0]),
            "sub": [1] * len(cssplits[0])}
    transposed_cssplits = [list(t) for t in zip(*cssplits)]
    for i, transposed_cssplit in enumerate(transposed_cssplits):
        for cs in transposed_cssplit:
            if cs.startswith("=") or cs == "N" or re.search(r"a|c|g|t|n", cs):
                continue
            if cs.startswith("+"):
                count["ins"][i] += len(cs.split("|"))
            elif cs.startswith("-"):
                count["del"][i] += 1
            elif cs.startswith("*"):
                count["sub"][i] += 1
    return count

def _remove_minor_indels(count_indels: dict[str, list[int]], coverage: int) -> dict[str, list[int]]:
    count_indels_removed = dict()
    threshold = coverage * 0.01
    for key, values in count_indels.items():
        values_removed = [v if v >= threshold else 1 for v in values]
        count_indels_removed.update({key: values_removed})
    return count_indels_removed


def _extract_anomaly_loci(indels_sample: dict, indels_control: dict) -> dict[str, set[int]]:
    anomaly_loci = dict()
    for key in indels_sample.keys():
        clf = LocalOutlierFactor(novelty=True)
        clf.fit(indels_control[key])
        pred = clf.predict(indels_sample[key])
        # はじめの9塩基はまず正常と判定されるはずなので、これをもとにして正常か異常かを判定する
        normal, abnormal = Counter(pred).keys()
        if Counter(pred[:9]).most_common()[0][0] == abnormal:
            normal, abnormal = abnormal, normal
        loci = {i for i, p in enumerate(pred) if p == abnormal}
        anomaly_loci.update({key: loci})
    return anomaly_loci

# def _extract_anomaly_loci(indels_sample: dict) -> dict[str, set[int]]:
#     anomaly_loci = dict()
#     for key, values in indels_sample.items():
#         clf = GaussianMixture(n_components = 2, random_state=0)
#         pred = clf.fit_predict(np.array(values).reshape(-1, 1))
#         # はじめの9塩基はまず正常と判定されるはずなので、これをもとにして正常か異常かを判定する
#         normal, abnormal = Counter(pred).keys()
#         if Counter(pred[:9]).most_common()[0][0] == abnormal:
#             normal, abnormal = abnormal, normal
#         loci = {i for i, p in enumerate(pred) if p == abnormal}
#         anomaly_loci.update({key: loci})
#     return anomaly_loci


# def _score_anomaly(count_sample: list[int | float]) -> list[float]:
#     x = np.array(count_sample)
#     average = np.mean(x)
#     sigma = np.std(x)
#     score = np.square((x - average) / sigma)
#     return score.tolist()


# def _extract_anomaly_loci_hotelling(score: list[float], threshold: float = 0.99) -> set[int]:
#     anomaly_loci = set()
#     thres = chi2.ppf(q=threshold, df=1)
#     for i, s in enumerate(score):
#         if s > thres:
#             anomaly_loci.add(i)
#     return anomaly_loci


def _split_kmer(indels: dict[str, list[int]], kmer: int = 10) -> dict[str, list[list[int]]]:
    results = defaultdict(list)
    center = kmer // 2
    for mut, value in indels.items():
        for i in range(len(value)):
            if center <= i <= len(value) - center:
                start = i - center
                if kmer % 2 == 0:
                    end = i + center
                else:
                    end = i + center + 1
                results[mut].append(value[start : end])
            else:
                results[mut].append([0]*kmer)
    return results


def _calc_distance(indels_sample: dict[str, list[list[int]]], indels_control: dict[str, list[list[int]]]) -> dict[str, list[float]]:
    results = defaultdict(list)
    for mut, value in indels_sample.items():
        for i, val in enumerate(value):
            dist = distance.euclidean(val, indels_control[mut][i])
            results[mut].append(dist)
    return results

def _extract(cssplits_sample, cssplits_control) -> dict[str, set[int]]:
    indels_sample = _count_indels(cssplits_sample)
    indels_control = _count_indels(cssplits_control)
    indels_sample = _remove_minor_indels(indels_sample, len(cssplits_sample))
    indels_control = _remove_minor_indels(indels_control, len(cssplits_control))
    # anomaly at single locus
    # mutation_locus_single = _extract_anomaly_loci(indels_sample)
    # mutation_locus_single = {}
    # loci_sample = _extract_anomaly_loci(indels_sample)
    # loci_control = _extract_anomaly_loci(indels_control)
    # for key in ["ins", "del", "sub"]:
    #     mutation_locus_single.update({key: loci_sample[key] - loci_control[key]})
    # Difference of anomaly within kmers
    x = _split_kmer(indels_sample, kmer = 10)
    y = _split_kmer(indels_control, kmer = 10)
    return _extract_anomaly_loci(x, y)
    # dists = _calc_distance(x, y)
    # loci_sample = _extract_anomaly_loci(dists)
    # mutation_locus_kmer = {}
    # for key in ["ins", "del", "sub"]:
    #     # score_sample = _score_anomaly(dists[key])
    #     mutation_locus_kmer.update({key: loci_sample[key]})
    # # Output
    # mutation_loci = {}
    # for key in ["ins", "del", "sub"]:
    #     mutation_loci.update({key: mutation_locus_single[key] & mutation_locus_kmer[key]})
    # return mutation_loci

###########################################################
# postprocesss
###########################################################


def _replace_errors_to_atmark(cssplits_sample: list[list[str]], mutation_loci: dict[str, set[int]]) -> list[list[str]]:
    results = []
    for cssplits in cssplits_sample:
        cssplits_replaced = []
        for i, cs in enumerate(cssplits):
            if cs.startswith("=") or cs == "N" or re.search(r"a|c|g|t|n", cs):
                cssplits_replaced.append(cs)
                continue
            if cs.startswith("+"):
                tag = "ins"
            elif cs.startswith("-"):
                tag = "del"
            else:
                tag = "sub"
            if i in mutation_loci[tag]:
                cssplits_replaced.append(cs)
            else:
                cssplits_replaced.append("@")
        results.append(cssplits_replaced)
    return results


def _replace_atmark(cssplits: list[list[str]], sequence: str) -> list[list[str]]:
    random.seed(1)
    cssplits_replaced = cssplits.copy()
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
            cs[0] = "N"
        if cs[-1] == "@":
            cs[-1] = "N"
    return cssplits_replaced

###############################################################################
# main
###############################################################################


def execute(TEMPDIR: Path, FASTA_ALLELES: dict[str, str], CONTROL_NAME: str, SAMPLE_NAME: str) -> None:
    for allele, sequence in FASTA_ALLELES.items():
        midsv_sample = midsv.read_jsonl((Path(TEMPDIR, "midsv", f"{SAMPLE_NAME}_{allele}.jsonl")))
        midsv_control = midsv.read_jsonl((Path(TEMPDIR, "midsv", f"{CONTROL_NAME}_{allele}.jsonl")))
        cssplits_sample = [cs["CSSPLIT"].split(",") for cs in midsv_sample]
        cssplits_control = [cs["CSSPLIT"].split(",") for cs in midsv_control]
        # Extract mutation loci
        mutation_loci = _extract(cssplits_sample, cssplits_control)
        # Correct sequence errors
        cssplits_sample_atmark = _replace_errors_to_atmark(cssplits_sample, mutation_loci)
        cssplits_control_atmark = _replace_errors_to_atmark(cssplits_control, mutation_loci)
        cssplits_sample_atmark_replaced = _replace_atmark(cssplits_sample_atmark, sequence)
        cssplits_control_atmark_replaced = _replace_atmark(cssplits_control_atmark, sequence)
        # Replace CSSPLIT
        cssplits_sample_corrected = [",".join(cs) for cs in cssplits_sample_atmark_replaced]
        cssplits_control_corrected = [",".join(cs) for cs in cssplits_control_atmark_replaced]
        for i, cssplits in enumerate(cssplits_sample_corrected):
            midsv_sample[i]["CSSPLIT"] = cssplits
        for i, cssplits in enumerate(cssplits_control_corrected):
            midsv_control[i]["CSSPLIT"] = cssplits
        midsv.write_jsonl(midsv_sample, Path(TEMPDIR, "midsv", f"{SAMPLE_NAME}_{allele}.jsonl"))
        midsv.write_jsonl(midsv_control, Path(TEMPDIR, "midsv", f"{CONTROL_NAME}_{allele}.jsonl"))

###############################################################################

# def _set_indexes(sequence: str):
#     sequence_length = len(sequence)
#     num_subset = sequence_length % 5
#     left_idx = 0
#     right_idx = sequence_length
#     if num_subset == 1:
#         left_idx += 1
#     elif num_subset == 2:
#         left_idx += 1
#         right_idx -= 1
#     elif num_subset == 3:
#         left_idx += 2
#         right_idx -= 1
#     elif num_subset == 4:
#         left_idx += 2
#         right_idx -= 2
#     return left_idx, right_idx


# def _count_5mer_indels(cssplits: list[list[str]], left_idx: int, right_idx: int) -> list[dict]:
#     transposed = [list(t) for t in zip(*cssplits)]
#     count_5mer = []
#     for i in range(left_idx, right_idx, 5):
#         count = {"ins": [1] * 5, "del": [1] * 5, "sub": [1] * 5}
#         cssplits_5mer = transposed[i : i + 5]
#         for j, cs in enumerate(cssplits_5mer):
#             counter = Counter(cs)
#             for key, cnt in counter.items():
#                 if key.startswith("=") or key == "N" or re.search(r"a|c|g|t|n", key):
#                     continue
#                 if key.startswith("+"):
#                     count["ins"][j] += cnt
#                 elif key.startswith("-"):
#                     count["del"][j] += cnt
#                 elif key.startswith("*"):
#                     count["sub"][j] += cnt
#         count_5mer.append(count)
#     return count_5mer


# def _remove_minor_indels(cssplits: list[list[str]], count_5mer: list[dict]) -> list[dict]:
#     coverage = len(cssplits)
#     count_5mer_filtered = []
#     for count in count_5mer:
#         dict_mutation = defaultdict(list)
#         for mutation in ["ins", "del", "sub"]:
#             if all(True if c < coverage * 0.01 else False for c in count[mutation]):
#                 dict_mutation[mutation] = [1] * 5
#             else:
#                 dict_mutation[mutation] = count[mutation]
#         count_5mer_filtered.append(dict_mutation)
#     return count_5mer_filtered


# def _extract_sequence_errors(count_5mer_sample, count_5mer_control):
#     sequence_errors = [set() for _ in range(len(count_5mer_sample))]
#     dists = defaultdict(list)
#     # Calculate Jensen-Shannon distance
#     for samp, cont in zip(count_5mer_sample, count_5mer_control):
#         for mutation in ["ins", "del", "sub"]:
#             s = samp[mutation]
#             c = cont[mutation]
#             dists[mutation].append(distance.jensenshannon(s, c))
#     # Discrimitate seq errors and real mutation using Hotelling's T-squared distribution
#     dists_all = np.array(list(dists.values())).flatten()
#     avg = np.average(dists_all[~np.isnan(dists_all)])
#     var = np.var(dists_all[~np.isnan(dists_all)])
#     threshold = 0.05
#     for mutation in ["ins", "del", "sub"]:
#         dists_subset = dists[mutation]
#         scores = [(xi - avg) ** 2 / var for xi in dists_subset]
#         thres = stats.chi2.interval(1 - threshold, 1)[1]
#         for i, score in enumerate(scores):
#             # 'nan' means the two distributions have too different, so it could be a real mutation
#             if np.isnan(score):
#                 continue
#             if score < thres:
#                 sequence_errors[i].add(mutation)
#     return sequence_errors


# def _replace_errors_to_atmark(cssplits_sample, sequence_errors, left_idx, right_idx):
#     cssplits_replaced = []
#     for samp in cssplits_sample:
#         samp_replaced = samp.copy()
#         for idx_error, idx_5mer in enumerate(range(left_idx, right_idx, 5)):
#             samp_5mer = samp[idx_5mer : idx_5mer + 5]
#             error = sequence_errors[idx_error]
#             if "ins" in error:
#                 samp_5mer = ["@" if cs.startswith("+") else cs for cs in samp_5mer]
#             if "del" in error:
#                 samp_5mer = ["@" if cs.startswith("-") else cs for cs in samp_5mer]
#             if "sub" in error:
#                 samp_5mer = ["@" if cs.startswith("*") else cs for cs in samp_5mer]
#             samp_replaced[idx_5mer : idx_5mer + 5] = samp_5mer
#         cssplits_replaced.append(samp_replaced)
#     return cssplits_replaced


# def _replace_atmark(cssplits: list[list[str]], sequence: str) -> list[list[str]]:
#     random.seed(1)
#     cssplits_replaced = cssplits.copy()
#     sequence_length = len(sequence)
#     for i in range(1, sequence_length - 1):
#         cssplits_atmark = defaultdict(str)
#         cssplits_sampling_key = defaultdict(list)
#         cssplits_sampling_all = []
#         flag_all_atmark = True
#         for idx, cssplit in enumerate(cssplits):
#             key = ",".join([cssplit[i - 1], cssplit[i + 1]])
#             if cssplit[i] == "@":
#                 cssplits_atmark[idx] = key
#             else:
#                 cssplits_sampling_key[key].append(cssplit[i])
#                 cssplits_sampling_all.append(cssplit[i])
#                 flag_all_atmark = False
#         for idx, key in cssplits_atmark.items():
#             if flag_all_atmark:
#                 cssplits_replaced[idx][i] = "N"
#             elif cssplits_sampling_key[key]:
#                 cssplits_replaced[idx][i] = random.choice(cssplits_sampling_key[key])
#             else:
#                 cssplits_replaced[idx][i] = random.choice(cssplits_sampling_all)
#     for cs in cssplits_replaced:
#         if cs[0] == "@":
#             cs[0] = "N"
#         if cs[-1] == "@":
#             cs[-1] = "N"
#     return cssplits_replaced


# ###############################################################################
# # main
# ###############################################################################


# def execute(TEMPDIR: Path, FASTA_ALLELES: dict[str, str], CONTROL_NAME: str, SAMPLE_NAME: str) -> None:
#     for allele, sequence in FASTA_ALLELES.items():
#         midsv_sample = midsv.read_jsonl((Path(TEMPDIR, "midsv", f"{SAMPLE_NAME}_{allele}.jsonl")))
#         midsv_control = midsv.read_jsonl((Path(TEMPDIR, "midsv", f"{CONTROL_NAME}_{allele}.jsonl")))
#         cssplits_sample = [cs["CSSPLIT"].split(",") for cs in midsv_sample]
#         cssplits_control = [cs["CSSPLIT"].split(",") for cs in midsv_control]
#         # Extract sequence errors
#         left_idx, right_idx = _set_indexes(sequence)
#         count_5mer_sample = _count_5mer_indels(cssplits_sample, left_idx, right_idx)
#         count_5mer_control = _count_5mer_indels(cssplits_control, left_idx, right_idx)
#         count_5mer_sample = _remove_minor_indels(cssplits_sample, count_5mer_sample)
#         count_5mer_control = _remove_minor_indels(cssplits_control, count_5mer_control)
#         sequence_errors = _extract_sequence_errors(count_5mer_sample, count_5mer_control)
#         # Correct sequence errors
#         cssplits_sample_error_replaced = _replace_errors_to_atmark(cssplits_sample, sequence_errors, left_idx, right_idx)
#         cssplits_control_error_replaced = _replace_errors_to_atmark(
#             cssplits_control, sequence_errors, left_idx, right_idx
#         )
#         cssplits_sample_atmark_replaced = _replace_atmark(cssplits_sample_error_replaced, sequence)
#         cssplits_control_atmark_replaced = _replace_atmark(cssplits_control_error_replaced, sequence)
#         # Replace CSSPLIT
#         cssplits_sample_corrected = [",".join(cs) for cs in cssplits_sample_atmark_replaced]
#         cssplits_control_corrected = [",".join(cs) for cs in cssplits_control_atmark_replaced]
#         for i, cssplits in enumerate(cssplits_sample_corrected):
#             midsv_sample[i]["CSSPLIT"] = cssplits
#         for i, cssplits in enumerate(cssplits_control_corrected):
#             midsv_control[i]["CSSPLIT"] = cssplits
#         midsv.write_jsonl(midsv_sample, Path(TEMPDIR, "midsv", f"{SAMPLE_NAME}_{allele}.jsonl"))
#         midsv.write_jsonl(midsv_control, Path(TEMPDIR, "midsv", f"{CONTROL_NAME}_{allele}.jsonl"))
