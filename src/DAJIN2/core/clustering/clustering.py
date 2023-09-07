from __future__ import annotations

import midsv
import random
from pathlib import Path
from itertools import groupby
from collections import defaultdict

from DAJIN2.core.clustering.score_handler import make_score, annotate_score
from DAJIN2.core.clustering.return_labels import return_labels
from DAJIN2.core.clustering.label_handler import relabel_with_consective_order

from DAJIN2.utils import io


# Constants
RANDOM_UPPER_LIMIT = 10**10
STRAND_BIAS_LOWER_LIMIT = 0.25
STRAND_BIAS_UPPER_LIMIT = 0.75


def is_strand_bias(path_control: Path) -> bool:
    count_strand = defaultdict(int)
    for m in midsv.read_jsonl(path_control):
        count_strand[m["STRAND"]] += 1

    total = count_strand["+"] + count_strand["-"]
    percentage_plus = count_strand["+"] / total if total else 0

    return not (STRAND_BIAS_LOWER_LIMIT < percentage_plus < STRAND_BIAS_UPPER_LIMIT)


def extract_labels(classif_sample, TEMPDIR, SAMPLE_NAME, CONTROL_NAME) -> list[dict[str]]:
    labels_all = []
    max_label = 0
    strand_bias = is_strand_bias(Path(TEMPDIR, CONTROL_NAME, "midsv", "control.json"))
    classif_sample.sort(key=lambda x: x["ALLELE"])
    for allele, group in groupby(classif_sample, key=lambda x: x["ALLELE"]):
        RANDOM_NUM = random.randint(0, 10**10)

        path_knockin_loci = Path(TEMPDIR, SAMPLE_NAME, "knockin_loci", f"{allele}.pickle")
        path_mutation_loci = Path(TEMPDIR, SAMPLE_NAME, "mutation_loci", f"{allele}.pickle")
        knockin_loci: set[int] = io.load_pickle(path_knockin_loci) if path_knockin_loci.exists() else set()
        mutation_loci: list[set[str]] = io.load_pickle(path_mutation_loci)
        if all(m == set() for m in mutation_loci):
            labels = [1] * len(classif_sample)
            labels_reorder = relabel_with_consective_order(labels, start=max_label)
            max_label = max(labels_reorder)
            labels_all.extend(labels_reorder)
            continue

        path_sample = Path(TEMPDIR, SAMPLE_NAME, "clustering", f"{allele}_{RANDOM_NUM}.json")
        path_control = Path(TEMPDIR, CONTROL_NAME, "midsv", f"{allele}.json")
        io.write_jsonl(data=group, file_path=path_sample)

        mutation_score: list[dict[str, float]] = make_score(path_sample, path_control, mutation_loci, knockin_loci)

        scores_sample = annotate_score(path_sample, mutation_score, mutation_loci)
        scores_control = annotate_score(path_control, mutation_score, mutation_loci, is_control=True)

        path_score_sample = Path(TEMPDIR, SAMPLE_NAME, "clustering", f"{allele}_score_{RANDOM_NUM}.json")
        path_score_control = Path(TEMPDIR, CONTROL_NAME, "clustering", f"{allele}_score_{RANDOM_NUM}.json")
        io.write_jsonl(data=scores_sample, file_path=path_score_sample)
        io.write_jsonl(data=scores_control, file_path=path_score_control)

        labels = return_labels(path_score_sample, path_score_control, path_sample, strand_bias)
        labels_reordered = relabel_with_consective_order(labels, start=max_label)

        max_label = max(labels_reordered)
        labels_all.extend(labels_reordered)
        # Remove temporary files
        path_sample.unlink()
        path_score_sample.unlink()
        path_score_control.unlink()

    return labels_all
