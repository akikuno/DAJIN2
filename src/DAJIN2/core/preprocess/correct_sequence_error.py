from __future__ import annotations

import random
import re
from collections import defaultdict
from pathlib import Path

import midsv
import numpy as np
from scipy import stats
from scipy.spatial import distance
from sklearn.neighbors import LocalOutlierFactor


# def _count_indels(cssplits: list[list[str]]) -> dict[str, list[int]]:
#     transposed_cssplits = [list(t) for t in zip(*cssplits)]
#     count = {"ins": [0] * len(transposed_cssplits),
#             "del": [0] * len(transposed_cssplits),
#             "sub": [0] * len(transposed_cssplits)}
#     for i, transposed_cssplit in enumerate(transposed_cssplits):
#         for cs in transposed_cssplit:
#             if cs.startswith("=") or cs == "N" or re.search(r"a|c|g|t|n", cs):
#                 continue
#             if cs.startswith("+"):
#                 # count["ins"][i] += len(cs.split("|"))
#                 count["ins"][i] += 1
#             elif cs.startswith("-"):
#                 count["del"][i] += 1
#             elif cs.startswith("*"):
#                 count["sub"][i] += 1
#     return count


# def _split_kmer(indels: dict[str, list[int]], kmer: int = 10) -> dict[str, list[list[int]]]:
#     results = defaultdict(list)
#     center = kmer // 2
#     for mut, value in indels.items():
#         for i in range(len(value)):
#             if center <= i <= len(value) - center:
#                 start = i - center
#                 if kmer % 2 == 0:
#                     end = i + center
#                 else:
#                     end = i + center + 1
#                 results[mut].append(value[start : end])
#             else:
#                 results[mut].append([0]*kmer)
#     return results


# def _extract_anomaly_loci(indels_kmer_sample: dict, indels_kmer_control: dict, coverage_sample: int, coverage_control: int) -> dict[str, set[int]]:
#     anomaly_loci = dict()
#     clf = LocalOutlierFactor(novelty=True, n_neighbors=5)
#     for key in indels_kmer_sample.keys():
#         loci = set()
#         values_control = np.array(indels_kmer_control[key]) / coverage_control
#         values_sample = np.array(indels_kmer_sample[key]) / coverage_sample
#         index = -1
#         for i, (value_control, value_sample) in enumerate(zip(values_control, values_sample)):
#             if i == index:
#                 continue
#             clf.fit(value_control.reshape(-1, 1))
#             pred = clf.predict(value_sample.reshape(-1, 1))
#             if pred[5] == -1:
#                 loci.add(i)
#             # If the next base is not -1, do not validate the next base because the next base is not an outlier.
#             if pred[6] == 1:
#                 index = i + 1
#         anomaly_loci.update({key: loci})
#     return anomaly_loci


# def _extract_dissimilar_loci(indels_kmer_sample: dict[str, list[list[int]]], indels_kmer_control: dict[str, list[list[int]]]) -> dict[str, set]:
#     results = dict()
#     for mut in indels_kmer_sample:
#         cossim = [distance.cosine(x, y) for x, y in zip(indels_kmer_sample[mut], indels_kmer_control[mut])]
#         pvalues = [stats.ttest_ind(x, y, equal_var=False)[1] for x, y in zip(indels_kmer_sample[mut], indels_kmer_control[mut])]
#         # if pvalue == nan, samples and controls are exactly same.
#         cossim_pval_false = [cossim if pvalue > 0.05 or np.isnan(pvalue) else 1 for cossim, pvalue in zip(cossim, pvalues)]
#         dissimilar_loci = {i for i, x in enumerate(cossim_pval_false) if x > 0.05}
#         results.update({mut: dissimilar_loci})
#     return results


# def _extract_mutation_loci(cssplits_sample, cssplits_control) -> dict[str, set[int]]:
#     indels_sample = _count_indels(cssplits_sample)
#     indels_control = _count_indels(cssplits_control)
#     indels_kmer_sample = _split_kmer(indels_sample, kmer = 10)
#     indels_kmer_control = _split_kmer(indels_control, kmer = 10)
#     anomaly_loci = _extract_anomaly_loci(indels_kmer_sample, indels_kmer_control, len(cssplits_sample), len(cssplits_control))
#     dissimilar_loci = _extract_dissimilar_loci(indels_kmer_sample, indels_kmer_control)
#     mutation_loci = dict()
#     for mut in anomaly_loci:
#         mutation_loci.update({mut: anomaly_loci[mut] & dissimilar_loci[mut]})
#     return mutation_loci

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


def correct_sequence_error(TEMPDIR: Path, FASTA_ALLELES: dict[str, str], CONTROL_NAME: str, SAMPLE_NAME: str, MUTATION_LOCI_ALLELES) -> None:
    for allele, sequence in FASTA_ALLELES.items():
        midsv_sample = midsv.read_jsonl((Path(TEMPDIR, "midsv", f"{SAMPLE_NAME}_{allele}.jsonl")))
        midsv_control = midsv.read_jsonl((Path(TEMPDIR, "midsv", f"{CONTROL_NAME}_{allele}.jsonl")))
        cssplits_sample = [cs["CSSPLIT"].split(",") for cs in midsv_sample]
        cssplits_control = [cs["CSSPLIT"].split(",") for cs in midsv_control]
        # Extract mutation loci
        mutation_loci = MUTATION_LOCI_ALLELES[allele]
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
