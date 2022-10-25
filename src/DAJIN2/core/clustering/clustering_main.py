from __future__ import annotations
import midsv
from pathlib import Path
from copy import deepcopy
from collections import defaultdict
from itertools import groupby
from .make_scores import make_scores
from .return_labels import return_labels
from .screen_diffloci import screen_different_loci


def extract_different_loci(TEMPDIR, classif_sample, MASKS_CONTROL, DICT_ALLELE, CONTROL_NAME):
    """
    Extract significantly different base loci between Sample and Control
    """
    dict_cssplit_control = defaultdict(list[dict])
    for allele in DICT_ALLELE.keys():
        path_control = Path(TEMPDIR, "midsv", f"{CONTROL_NAME}_{allele}.jsonl")
        cssplit_control = [cs["CSSPLIT"] for cs in midsv.read_jsonl(path_control)]
        dict_cssplit_control[allele] = cssplit_control
    classif_sample.sort(key=lambda x: (x["ALLELE"], x["SV"]))
    DIFFLOCI = defaultdict(list[dict])
    REPETITIVE_DELLOCI = defaultdict(list[dict])
    for (allele, sv), group in groupby(classif_sample, key=lambda x: (x["ALLELE"], x["SV"])):
        cssplit_sample = [record["CSSPLIT"] for record in group]
        cssplit_control = dict_cssplit_control[allele]
        sequence = DICT_ALLELE[allele]
        masks_control = MASKS_CONTROL[allele]
        diffloci, repetitive_delloci = screen_different_loci(
            cssplit_sample, cssplit_control, sequence, masks_control, alpha=0.01, threshold=0.05
        )
        DIFFLOCI[f"{allele}-{sv}"] = diffloci
        REPETITIVE_DELLOCI[f"{allele}-{sv}"] = repetitive_delloci
    return DIFFLOCI, REPETITIVE_DELLOCI


def add_labels(classif_sample, diffloci_by_alleles):
    labels = []
    label_start = 1
    classif_sample.sort(key=lambda x: (x["ALLELE"], x["SV"]))
    for (allele, sv), group in groupby(classif_sample, key=lambda x: (x["ALLELE"], x["SV"])):
        key = f"{allele}-{sv}"
        cssplit_sample = [g["CSSPLIT"] for g in group]
        diffloci = diffloci_by_alleles[key]
        scores = []
        if diffloci is not None:
            scores = make_scores(cssplit_sample, diffloci)
        if any(scores):
            labels += [label + label_start for label in return_labels(scores).tolist()]
        else:
            labels += [label_start] * len(cssplit_sample)
        label_start = len(set(labels)) + 1
    clust_sample = deepcopy(classif_sample)
    for clust, label in zip(clust_sample, labels):
        clust["LABEL"] = label
    return clust_sample


def add_readnum(clust_sample):
    readnum = defaultdict(int)
    for cs in clust_sample:
        readnum[cs["LABEL"]] += 1
    for cs in clust_sample:
        cs["READNUM"] = readnum[cs["LABEL"]]
    return clust_sample


def add_percent(clust_sample):
    n_sample = len(clust_sample)
    percent = defaultdict(int)
    for cs in clust_sample:
        percent[cs["LABEL"]] += 1 / n_sample
    percent = {key: round(val * 100, 3) for key, val in percent.items()}
    for cs in clust_sample:
        cs["PERCENT"] = percent[cs["LABEL"]]
    return clust_sample


def update_labels(clust_sample):
    """
    Allocate new labels according to the ranking by PERCENT
    """
    clust_sample.sort(key=lambda x: (-x["PERCENT"], x["LABEL"]))
    new_label = 1
    prev_label = clust_sample[0]["LABEL"]
    for cs in clust_sample:
        if prev_label != cs["LABEL"]:
            new_label += 1
        prev_label = cs["LABEL"]
        cs["LABEL"] = new_label
    return clust_sample
