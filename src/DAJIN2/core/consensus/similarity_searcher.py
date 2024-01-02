from __future__ import annotations

from pathlib import Path
from collections import defaultdict

import numpy as np

from sklearn.neighbors import LocalOutlierFactor

from DAJIN2.utils import io

"""
Sequence errors (such as large deletions) present in the control significantly impair the accuracy of mutation_loci. Therefore, only reads similar to alleles after clustering are desired for analysis.
Additionally, the mutation regions in the alleles after clustering should be identified and masked.
"""


def onehot_by_mutations(midsv_sample: list[dict]) -> dict[str, np.ndarray]:
    mut_onehot = defaultdict(list)
    for c in midsv_sample:
        cssplits = c["CSSPLIT"].split(",")
        for mut in {"+", "-", "*"}:
            onehot = [1 if cs.startswith(mut) else 0 for cs in cssplits]
            mut_onehot[mut].append(onehot)
    return {mut: np.array(value) for mut, value in mut_onehot.items()}


def calculate_percentage(
    mut_onehot_sample: dict[str, np.ndarray], coverage_match: np.ndarray[int]
) -> dict[str, np.ndarray]:
    mut_percentage = dict()
    for mut, onehot in mut_onehot_sample.items():
        x = np.sum(onehot, axis=0) / coverage_match
        mut_percentage[mut] = np.where(np.isnan(x), 0, x)
    return mut_percentage


def get_values_to_mask(mut_percentage_sample: dict[str, np.ndarray], threshold=0.5) -> dict[str, np.ndarray[float]]:
    mask = dict()
    for mut, percentage in mut_percentage_sample.items():
        mask[mut] = np.where(percentage > threshold, 0, percentage)
    return mask


def apply_mask(mut_onehot: dict[str, np.ndarray], mask_sample: dict[str, np.ndarray[float]]):
    mut_onehot_masked = dict()
    for mut, onehot in mut_onehot.items():
        mut_onehot_masked[mut] = onehot * mask_sample[mut]
    return mut_onehot_masked


def identify_normal_reads(
    mut_onehot_sample_masked: dict[str, np.ndarray], mut_onehot_control_masked: dict[str, np.ndarray]
) -> list[bool]:
    mutation_comparisons = dict()
    for mut in {"+", "-", "*"}:
        values_sample = mut_onehot_sample_masked[mut]
        values_control = mut_onehot_control_masked[mut]

        values_sum_sample = values_sample.sum(axis=1).reshape(-1, 1)
        values_sum_control = values_control.sum(axis=1).reshape(-1, 1)

        clf = LocalOutlierFactor(novelty=True, n_neighbors=len(values_sum_sample)).fit(values_sum_sample)
        labels = clf.predict(values_sum_control)
        mutation_comparisons[mut] = np.where(labels == 1, True, False)

    return (mutation_comparisons["+"] * mutation_comparisons["-"] * mutation_comparisons["*"]).tolist()


###########################################################
# main
###########################################################


def filter_control(path_midsv_control: Path, path_midsv_sample: Path) -> list[bool]:
    """
    find similar control reads compared to sample reads
    """
    cssplits = (m["CSSPLIT"].split(",") for m in io.read_jsonl(path_midsv_sample))
    coverage_match = np.array([sum(1 for cs in cssplit if cs.startswith("=")) for cssplit in zip(*cssplits)])
    mut_onehot_sample = onehot_by_mutations(io.read_jsonl(path_midsv_sample))
    mut_onehot_control = onehot_by_mutations(io.read_jsonl(path_midsv_control))

    mut_percentage_sample = calculate_percentage(mut_onehot_sample, coverage_match)
    values_mask = get_values_to_mask(mut_percentage_sample)

    mut_onehot_sample_masked = apply_mask(mut_onehot_sample, values_mask)
    mut_onehot_control_masked = apply_mask(mut_onehot_control, values_mask)

    return identify_normal_reads(mut_onehot_sample_masked, mut_onehot_control_masked)


def cache_selected_control_by_similarity(path_midsv_control: Path, path_midsv_sample: Path, path_output: Path) -> None:
    normal_reads_flags = filter_control(path_midsv_control, path_midsv_sample)
    midsv_control = io.read_jsonl(path_midsv_control)
    midsv_filtered = (m for m, flag in zip(midsv_control, normal_reads_flags) if flag is True)

    allele, label, *_ = Path(path_midsv_sample).stem.split("_")
    io.write_jsonl(midsv_filtered, Path(path_output, f"{allele}_{label}_control.jsonl"))
