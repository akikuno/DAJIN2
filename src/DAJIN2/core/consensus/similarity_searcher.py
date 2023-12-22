from __future__ import annotations

from pathlib import Path
from collections import defaultdict

import numpy as np
from sklearn.cluster import MiniBatchKMeans

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


def calc_percentage(mut_onehot: dict[str, np.ndarray], coverage: int) -> dict[str, np.ndarray]:
    mut_percentage = dict()
    for mut, onehot in mut_onehot.items():
        mut_percentage[mut] = np.sum(onehot, axis=0) / coverage
    return mut_percentage


def get_index_to_mask(mut_percentage_sample: dict[str, np.ndarray], threshold=0.5) -> dict[str, np.ndarray[bool]]:
    mask = dict()
    for mut, percentage in mut_percentage_sample.items():
        mask[mut] = np.where(percentage > threshold, 0, percentage)
    return mask


def apply_mask(mut_onehot: dict[str, np.ndarray], mask_sample: dict[str, np.ndarray[bool]]):
    mut_onehot_masked = dict()
    for mut, onehot in mut_onehot.items():
        mut_onehot_masked[mut] = onehot * mask_sample[mut]
    return mut_onehot_masked


def calc_similarity(mut_onehot_sample_masked: dict[str, np.ndarray], mut_onehot_control_masked: dict[str, np.ndarray]):
    """The label compares the number of mutations in the sample and control.
    The number of mutations in the sample is classified into two categories, and a threshold for being considered normal is determined, which is then applied to the control.

    On the other hand, the distance compares the Euclidean distance between the sample (average value of mutations at each index) and the control.
    The smaller the Euclidean distance, the more similar they are.
    """
    similarity = dict()
    for mut in {"+", "-", "*"}:
        values_sample = mut_onehot_sample_masked[mut]
        values_control = mut_onehot_control_masked[mut]

        kmeans = MiniBatchKMeans(n_clusters=2, random_state=0).fit(values_sample.sum(axis=1).reshape(-1, 1))
        threshold = kmeans.cluster_centers_.mean()
        label = np.where(values_control.sum(axis=1) <= threshold, 1, 0)

        euclidean_distance = np.linalg.norm(values_control - values_sample.mean(axis=0), axis=1)

        similarity[mut] = {"label": label, "distance": euclidean_distance}

    return similarity


def is_similar(similarity, coverage_sample) -> list[bool]:
    x = similarity["+"]["label"] * similarity["-"]["label"] * similarity["*"]["label"]
    if np.sum(x) >= coverage_sample:
        return np.where(x == 1, True, False).tolist()
    else:
        # If there are few similar reads, extract reads with a smaller distance
        x = similarity["+"]["distance"] + similarity["-"]["distance"] + similarity["*"]["distance"]
        n = min(coverage_sample, len(x))
        threshold = np.sort(x)[n - 1]
        return np.where(x < threshold, True, False).tolist()


###########################################################
# main
###########################################################


def filter_control(path_midsv_control: Path, path_midsv_sample: Path) -> list[bool]:
    """
    find similar control reads compared to sample reads
    """
    coverage_sample = io.count_newlines(path_midsv_sample)
    mut_onehot_sample = onehot_by_mutations(io.read_jsonl(path_midsv_sample))
    mut_onehot_control = onehot_by_mutations(io.read_jsonl(path_midsv_control))

    mut_percentage_sample = calc_percentage(mut_onehot_sample, coverage_sample)
    mask_sample = get_index_to_mask(mut_percentage_sample)

    mut_onehot_sample_masked = apply_mask(mut_onehot_sample, mask_sample)
    mut_onehot_control_masked = apply_mask(mut_onehot_control, mask_sample)

    similarity = calc_similarity(mut_onehot_sample_masked, mut_onehot_control_masked)

    return is_similar(similarity, coverage_sample)


def cache_selected_control_by_similarity(path_midsv_control: Path, path_midsv_sample: Path, path_output: Path) -> None:
    bools = filter_control(path_midsv_control, path_midsv_sample)
    midsv_control = io.read_jsonl(path_midsv_control)
    midsv_filtered = (m for m, b in zip(midsv_control, bools) if b is True)

    allele, label, *_ = Path(path_midsv_sample).stem.split("_")
    io.write_jsonl(midsv_filtered, Path(path_output, f"{allele}_{label}_control.jsonl"))
