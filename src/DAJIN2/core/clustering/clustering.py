from __future__ import annotations
from pathlib import Path
from itertools import groupby
import midsv
from copy import deepcopy
from collections import defaultdict
from src.DAJIN2.core.clustering.correct_cssplits import correct_cssplits
from src.DAJIN2.core.clustering.make_score import make_score
from src.DAJIN2.core.clustering.annotate_score import annotate_score
from src.DAJIN2.core.clustering.merge_clusters import merge_clusters
from src.DAJIN2.core.clustering.reorder_labels import reorder_labels
from src.DAJIN2.core.clustering.return_labels import return_labels

# find_knockin_loci


def add_labels(classif_sample, classif_control, FASTA_ALLELES: dict, THREADS: int = 1) -> list[dict[str]]:
    labels_all = []
    max_label = 0
    # KNOCKIN_LOCI = find_knockin_loci(TEMPDIR, FASTA_ALLELES, CONTROL_NAME)
    classif_sample.sort(key=lambda x: (x["ALLELE"], x["SV"]))
    for (allele, _), group in groupby(classif_sample, key=lambda x: (x["ALLELE"], x["SV"])):
        # sequence = FASTA_ALLELES[allele]
        # knockin_loci = KNOCKIN_LOCI[allele]
        # Control
        cssplits_control = [
            cs["CSSPLIT"].split(",") for cs in classif_control if cs["ALLELE"] == allele and cs["SV"] == False
        ]
        # Sample
        cssplits_sample = [cs["CSSPLIT"].split(",") for cs in group]
        mutation_score = make_score(cssplits_control, cssplits_sample)
        scores_control = annotate_score(cssplits_control, mutation_score)
        scores_sample = annotate_score(cssplits_sample, mutation_score)
        scores = scores_sample + scores_control[:1000]
        labels = return_labels(scores, THREADS)
        labels_control = labels[len(scores_sample) :]
        labels_sample = labels[: len(scores_sample)]
        labels_merged = merge_clusters(labels_control, labels_sample)
        labels_reorder = reorder_labels(labels_merged, start=max_label)
        max_label = max(labels_reorder)
        labels_all.extend(labels_reorder)
    clust_sample = deepcopy(classif_sample)
    for clust, label in zip(clust_sample, labels_all):
        clust["LABEL"] = label
    return clust_sample


def add_readnum(clust_sample):
    clust_result = deepcopy(clust_sample)
    readnum = defaultdict(int)
    for cs in clust_result:
        readnum[cs["LABEL"]] += 1
    for cs in clust_result:
        cs["READNUM"] = readnum[cs["LABEL"]]
    return clust_result


def add_percent(clust_sample):
    clust_result = deepcopy(clust_sample)
    n_sample = len(clust_result)
    percent = defaultdict(int)
    for cs in clust_result:
        percent[cs["LABEL"]] += 1 / n_sample
    percent = {key: round(val * 100, 3) for key, val in percent.items()}
    for cs in clust_result:
        cs["PERCENT"] = percent[cs["LABEL"]]
    return clust_result


def update_labels(clust_sample):
    """
    Allocate new labels according to the ranking by PERCENT
    """
    clust_result = deepcopy(clust_sample)
    clust_result.sort(key=lambda x: (-x["PERCENT"], x["LABEL"]))
    new_label = 1
    prev_label = clust_result[0]["LABEL"]
    for cs in clust_result:
        if prev_label != cs["LABEL"]:
            new_label += 1
        prev_label = cs["LABEL"]
        cs["LABEL"] = new_label
    return clust_result
