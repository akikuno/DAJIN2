from __future__ import annotations

from DAJIN2.utils import io, config

config.set_warnings()

import numpy as np

from pathlib import Path
from collections import Counter
from collections import defaultdict

from sklearn.mixture import GaussianMixture

import random
from itertools import groupby

from DAJIN2.core.clustering.score_handler import make_score, annotate_score
from DAJIN2.core.clustering.label_handler import relabel_with_consective_order
from DAJIN2.core.clustering.clustering import reduce_dimension, optimize_labels, get_labels_removed_strand_bias


# Constants
RANDOM_UPPER_LIMIT = 10**10
STRAND_BIAS_LOWER_LIMIT = 0.25
STRAND_BIAS_UPPER_LIMIT = 0.75


def is_strand_bias(path_control: Path) -> bool:
    count_strand = defaultdict(int)
    for m in io.read_jsonl(path_control):
        count_strand[m["STRAND"]] += 1

    total = count_strand["+"] + count_strand["-"]
    percentage_plus = count_strand["+"] / total if total else 0

    return not (STRAND_BIAS_LOWER_LIMIT < percentage_plus < STRAND_BIAS_UPPER_LIMIT)


def get_label_most_common(labels: list[int]) -> int:
    return Counter(labels).most_common()[0][0]


def subset_scores(labels: list[int], scores: list[int], label_most: int, length: int = 1000) -> list[int]:
    subset = [score for label, score in zip(labels, scores) if label == label_most]
    return subset[:length]


def return_labels(
    path_score_sample: Path | str, path_score_control: Path | str, path_sample: Path | str, strand_bias: bool
) -> list[int]:
    np.random.seed(seed=1)
    X_control = reduce_dimension([], io.read_jsonl(path_score_control))
    # subset to 1000 reads of controls in the most common cluster to remove outliers and reduce computation time
    labels_control = GaussianMixture(n_components=2, random_state=1).fit_predict(X_control)
    label_most_common = get_label_most_common(labels_control)
    scores_control_subset = subset_scores(labels_control, io.read_jsonl(path_score_control), label_most_common, 1000)
    X = reduce_dimension(io.read_jsonl(path_score_sample), scores_control_subset)
    coverage_sample = io.count_newlines(path_score_sample)
    coverage_control = len(scores_control_subset)
    labels = optimize_labels(X, coverage_sample, coverage_control)
    # correct clusters with strand bias
    if strand_bias is False:
        labels = get_labels_removed_strand_bias(path_sample, path_score_sample, labels)
    return labels


def extract_labels(classif_sample, TEMPDIR, SAMPLE_NAME, CONTROL_NAME) -> list[dict[str]]:
    labels_all = []
    max_label = 0
    strand_bias = is_strand_bias(Path(TEMPDIR, CONTROL_NAME, "midsv", "control.json"))
    classif_sample.sort(key=lambda x: x["ALLELE"])
    for allele, group in groupby(classif_sample, key=lambda x: x["ALLELE"]):
        path_mutation_loci = Path(TEMPDIR, SAMPLE_NAME, "mutation_loci", f"{allele}.pickle")
        mutation_loci: list[set[str]] = io.load_pickle(path_mutation_loci)
        if all(m == set() for m in mutation_loci):
            max_label += 1
            labels_all.extend([max_label] * len(classif_sample))
            continue

        path_knockin_loci = Path(TEMPDIR, SAMPLE_NAME, "knockin_loci", f"{allele}.pickle")
        knockin_loci: set[int] = io.load_pickle(path_knockin_loci) if path_knockin_loci.exists() else set()

        RANDOM_INT = random.randint(0, RANDOM_UPPER_LIMIT)
        path_sample = Path(TEMPDIR, SAMPLE_NAME, "clustering", f"{allele}_{RANDOM_INT}.json")
        path_control = Path(TEMPDIR, CONTROL_NAME, "midsv", f"{allele}.json")
        io.write_jsonl(data=group, file_path=path_sample)

        """Prepare and write clustering data to temporary files."""
        mutation_score: list[dict[str, float]] = make_score(path_sample, path_control, mutation_loci, knockin_loci)

        scores_sample = annotate_score(path_sample, mutation_score, mutation_loci)
        scores_control = annotate_score(path_control, mutation_score, mutation_loci, is_control=True)

        path_score_sample = Path(TEMPDIR, SAMPLE_NAME, "clustering", f"{allele}_score_{RANDOM_INT}.json")
        path_score_control = Path(TEMPDIR, CONTROL_NAME, "clustering", f"{allele}_score_{RANDOM_INT}.json")
        io.write_jsonl(data=scores_sample, file_path=path_score_sample)
        io.write_jsonl(data=scores_control, file_path=path_score_control)

        """Extract labels."""
        labels = return_labels(path_score_sample, path_score_control, path_sample, strand_bias)
        labels_reordered = relabel_with_consective_order(labels, start=max_label)

        max_label = max(labels_reordered)
        labels_all.extend(labels_reordered)

        """Remove temporary files."""
        path_sample.unlink()
        path_score_sample.unlink()
        path_score_control.unlink()

    return labels_all
