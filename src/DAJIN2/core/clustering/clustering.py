from __future__ import annotations

from collections import defaultdict
from itertools import groupby
from pathlib import Path

import midsv

from DAJIN2.core.clustering.make_score import make_score
from DAJIN2.core.clustering.preprocess import (compress_insertion,
                                               replace_both_ends_n)
from DAJIN2.core.clustering.return_labels import return_labels
from DAJIN2.core.preprocess.correct_knockin import extract_knockin_loci


def extract_cssplits_in_mutation(cssplits_sample: list[list], mutation_loci: set) -> list[list]:
    cssplits_mutation = []
    for cssplits in cssplits_sample:
        cs_mutation = []
        for i, cs in enumerate(cssplits):
            if i in mutation_loci:
                cs_mutation.append(cs)
        cssplits_mutation.append(cs_mutation)
    return cssplits_mutation


def annotate_score(cssplits: list[list], mutation_score: list[dict]):
    scores = []
    for cssplit in cssplits:
        score = [0]
        for i in range(1, len(cssplit) - 1):
            if not mutation_score[i]:
                score.append(0)
                continue
            kmer = ",".join([cssplit[i - 1], cssplit[i], cssplit[i + 1]])
            score.append(mutation_score[i].get(kmer, 0))
        scores.append(score + [0])
    return scores


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


def add_labels(classif_sample, TEMPDIR, CONTROL_NAME, MUTATION_LOCI, THREADS: int = 1) -> list[dict[str]]:
    paths_midsv = list(Path(TEMPDIR, "midsv").glob(f"{CONTROL_NAME}_splice_*"))
    cssplits_control_by_alleles = defaultdict(list)
    for path_midsv in paths_midsv:
        midsv_control = midsv.read_jsonl(path_midsv)
        allele = path_midsv.stem.split("_")[-1]
        cssplits = [cs["CSSPLIT"].split(",") for cs in midsv_control]
        cssplits_control_by_alleles[allele] = cssplits
    knockin_alleles = extract_knockin_loci(TEMPDIR)
    labels_all = []
    max_label = 0
    classif_sample.sort(key=lambda x: x["ALLELE"])
    for allele, group in groupby(classif_sample, key=lambda x: x["ALLELE"]):
        mutation_loci: set = MUTATION_LOCI[allele]
        cssplits_control = cssplits_control_by_alleles[allele]
        cssplits_sample = [cs["CSSPLIT"].split(",") for cs in group]
        # cssplits_control = replace_both_ends_n(cssplits_control)
        # cssplits_sample = replace_both_ends_n(cssplits_sample)
        cssplits_control = extract_cssplits_in_mutation(cssplits_control, mutation_loci)
        cssplits_sample = extract_cssplits_in_mutation(cssplits_sample, mutation_loci)
        cssplits_control = compress_insertion(cssplits_control)
        cssplits_sample = compress_insertion(cssplits_sample)
        mutation_score = make_score(cssplits_control, cssplits_sample, knockin_alleles[allele])
        scores_control = annotate_score(cssplits_control, mutation_score)
        scores_sample = annotate_score(cssplits_sample, mutation_score)
        labels = return_labels(scores_sample, scores_control)
        labels_reorder = reorder_labels(labels, start=max_label)
        max_label = max(labels_reorder)
        labels_all.extend(labels_reorder)
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
