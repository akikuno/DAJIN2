from __future__ import annotations

import pickle
import midsv
import random
from pathlib import Path
from itertools import groupby
from typing import Generator
from collections import defaultdict

from DAJIN2.core.clustering.make_kmer import generate_mutation_kmers
from DAJIN2.core.clustering.make_score import make_score
from DAJIN2.core.clustering.return_labels import return_labels
from DAJIN2.utils import io


def annotate_score(path_sample, mutation_score, mutation_loci, is_control=False) -> Generator[list[float]]:
    for cssplit_kmer in generate_mutation_kmers(path_sample, mutation_loci):
        score = [0 for _ in range(len(cssplit_kmer))]
        for i, (cs_kmer, mut_score) in enumerate(zip(cssplit_kmer, mutation_score)):
            if mut_score == {}:
                continue
            # Mutation sites are not considered in controls
            # because they should be sample-specific.
            if is_control and cs_kmer.split(",")[1][0] in mutation_loci[i]:
                continue
            if cs_kmer in mut_score:
                score[i] = mut_score[cs_kmer]
        yield score


def reorder_labels(labels: list[int], start: int = 0) -> list[int]:
    labels_ordered = labels.copy()
    num = start
    d = defaultdict(int)
    for i, l in enumerate(labels_ordered):
        if not d[l]:
            num += 1
            d[l] = num
        labels_ordered[i] = d[l]
    return labels_ordered


###########################################################
# main
###########################################################


def is_strand_bias(path_control) -> bool:
    count_strand = defaultdict(int)
    for m in midsv.read_jsonl(path_control):
        count_strand[m["STRAND"]] += 1
    percentage_plus = count_strand["+"] / (count_strand["+"] + count_strand["-"])
    if 0.25 < percentage_plus < 0.75:
        return False
    else:
        return True


def add_labels(classif_sample, TEMPDIR, SAMPLE_NAME, CONTROL_NAME) -> list[dict[str]]:
    labels_all = []
    max_label = 0
    strand_bias = is_strand_bias(Path(TEMPDIR, CONTROL_NAME, "midsv", "control.json"))
    classif_sample.sort(key=lambda x: x["ALLELE"])
    for allele, group in groupby(classif_sample, key=lambda x: x["ALLELE"]):
        RANDOM_NUM = random.randint(0, 10**10)
        path_knockin = Path(TEMPDIR, SAMPLE_NAME, "knockin_loci", f"{allele}.pickle")
        if path_knockin.exists():
            with open(path_knockin, "rb") as p:
                knockin_loci = pickle.load(p)
        else:
            knockin_loci = set()
        with open(Path(TEMPDIR, SAMPLE_NAME, "mutation_loci", f"{allele}.pickle"), "rb") as p:
            mutation_loci: dict[str, set[str]] = pickle.load(p)
        if all(m == set() for m in mutation_loci):
            labels = [1] * len(classif_sample)
            labels_reorder = reorder_labels(labels, start=max_label)
            max_label = max(labels_reorder)
            labels_all.extend(labels_reorder)
            continue
        path_sample = Path(TEMPDIR, SAMPLE_NAME, "clustering", f"{allele}_{RANDOM_NUM}.json")
        path_control = Path(TEMPDIR, CONTROL_NAME, "midsv", f"{allele}.json")
        io.write_jsonl(data=group, path=path_sample)
        mutation_score: list[dict[str, float]] = make_score(path_sample, path_control, mutation_loci, knockin_loci)
        scores_sample = annotate_score(path_sample, mutation_score, mutation_loci)
        scores_control = annotate_score(path_control, mutation_score, mutation_loci, is_control=True)
        path_score_sample = Path(TEMPDIR, SAMPLE_NAME, "clustering", f"{allele}_score_{RANDOM_NUM}.json")
        path_score_control = Path(TEMPDIR, CONTROL_NAME, "clustering", f"{allele}_score_{RANDOM_NUM}.json")
        io.write_jsonl(data=scores_sample, path=path_score_sample)
        io.write_jsonl(data=scores_control, path=path_score_control)
        labels = return_labels(path_score_sample, path_score_control, path_sample, strand_bias)
        labels_reorder = reorder_labels(labels, start=max_label)
        max_label = max(labels_reorder)
        labels_all.extend(labels_reorder)
        # Remove temporary files
        path_sample.unlink()
        path_score_sample.unlink()
        path_score_control.unlink()
    clust_sample = classif_sample.copy()
    for clust, label in zip(clust_sample, labels_all):
        clust["LABEL"] = label
    return clust_sample


def add_readnum(clust_sample: list[dict]) -> list[dict]:
    clust_result = clust_sample.copy()
    readnum = defaultdict(int)
    for cs in clust_result:
        readnum[cs["LABEL"]] += 1
    for cs in clust_result:
        cs["READNUM"] = readnum[cs["LABEL"]]
    return clust_result


def add_percent(clust_sample: list[dict]) -> list[dict]:
    clust_result = clust_sample.copy()
    n_sample = len(clust_result)
    percent = defaultdict(int)
    for cs in clust_result:
        percent[cs["LABEL"]] += 1 / n_sample
    percent = {key: round(val * 100, 3) for key, val in percent.items()}
    for cs in clust_result:
        cs["PERCENT"] = percent[cs["LABEL"]]
    return clust_result


def update_labels(clust_sample: list[dict]) -> list[dict]:
    """
    Allocate new labels according to the ranking by PERCENT
    """
    clust_result = clust_sample.copy()
    clust_result.sort(key=lambda x: (-x["PERCENT"], x["LABEL"]))
    new_label = 1
    prev_label = clust_result[0]["LABEL"]
    for cs in clust_result:
        if prev_label != cs["LABEL"]:
            new_label += 1
        prev_label = cs["LABEL"]
        cs["LABEL"] = new_label
    return clust_result
