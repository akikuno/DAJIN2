from __future__ import annotations

from collections import defaultdict
from itertools import groupby
from pathlib import Path

import midsv

from DAJIN2.core.clustering.make_score import make_score
from DAJIN2.core.clustering.return_labels import return_labels


def _compress_insertion(cssplits: list[list[str]]) -> list[dict[str, int]]:
    """Insertion will be subdivided by sequence error in the its sequence, so it is compressed as a '+I' to eliminate mutations.
    #TODO ただ、これでは、insertion配列の中に真のmutationがある場合に、そのmutationを抽出できないので、**insertion配列の中にmutationがある場合は、insertion配列をそのまま残す**必要がある。
    """
    cssplits_abstracted = []
    for cssplit in cssplits:
        for i, cs in enumerate(cssplit):
            if cs.startswith("+"):
                cssplit[i] = "+I" + cs.split("|")[-1]
        cssplits_abstracted.append(cssplit)
    return cssplits_abstracted


def _transpose_mutation_loci(mutation_loci, cssplits_sample):
    mutation_loci_transposed = [set() for _ in range(len(cssplits_sample[0]))]
    for mut, values in mutation_loci.items():
        for i, loci in enumerate(mutation_loci_transposed):
            if i in values:
                loci.add(mut)
    return mutation_loci_transposed


def _extract_cssplits_in_mutation_by_3mer(cssplits: list[list[str]], mutation_loci_transposed: list[set[str]]) -> list[list[str]]:
    cssplits_mutation = []
    for cssplit in cssplits:
        cs_mutation = ["N,N,N"]
        for i in range(1, len(cssplit) - 1):
            if mutation_loci_transposed[i] == set():
                cs_mutation.append("N,N,N")
                continue
            mutation = ""
            if cssplit[i].startswith("+"):
                mutation = "ins"
            elif cssplit[i].startswith("-"):
                mutation = "del"
            elif cssplit[i].startswith("*"):
                mutation = "sub"
            if mutation in mutation_loci_transposed[i]:
                kmer = ",".join([cssplit[i - 1], cssplit[i], cssplit[i + 1]])
                cs_mutation.append(kmer)
            else:
                cs_mutation.append("N,N,N")
        cs_mutation.append("N,N,N")
        cssplits_mutation.append(cs_mutation)
    return cssplits_mutation


def _annotate_score(cssplits: list[list[str]], mutation_score: list[dict[str:float]]) -> list[list[float]]:
    scores = []
    for cssplit in cssplits:
        score = [0 for _ in range(len(cssplit))]
        for i, (cs, mutscore) in enumerate(zip(cssplit, mutation_score)):
            if mutscore == {}:
                continue
            if cs in mutscore:
                score[i] = mutscore[cs]
        scores.append(score)
    return scores


def _reorder_labels(labels: list[int], start: int = 0) -> list[int]:
    labels_ordered = labels.copy()
    num = start
    d = defaultdict(int)
    for i, l in enumerate(labels_ordered):
        if not d[l]:
            num += 1
            d[l] = num
        labels_ordered[i] = d[l]
    return labels_ordered


def add_labels(classif_sample, TEMPDIR, CONTROL_NAME, MUTATION_LOCI_ALLELES, KNOCKIN_LOCI_ALLELES, THREADS: int = 1) -> list[dict[str]]:
    labels_all = []
    max_label = 0
    classif_sample.sort(key=lambda x: x["ALLELE"])
    for allele, group in groupby(classif_sample, key=lambda x: x["ALLELE"]):
        mutation_loci: dict[str, set[int]] = MUTATION_LOCI_ALLELES[allele]
        if sum(len(v) for v in mutation_loci.values()) == 0:
            labels = [1] * len(classif_sample)
            labels_reorder = _reorder_labels(labels, start=max_label)
            max_label = max(labels_reorder)
            labels_all.extend(labels_reorder)
            continue
        knockin_loci = set()
        if allele in KNOCKIN_LOCI_ALLELES:
            knockin_loci = KNOCKIN_LOCI_ALLELES[allele]
        cssplits_sample = [cs["CSSPLIT"].split(",") for cs in group]
        midsv_control = midsv.read_jsonl((Path(TEMPDIR, "midsv", f"{CONTROL_NAME}_{allele}.jsonl")))
        cssplits_control = [cs["CSSPLIT"].split(",") for cs in midsv_control]
        cssplits_sample = _compress_insertion(cssplits_sample)
        cssplits_control = _compress_insertion(cssplits_control)
        mutation_loci_transposed = _transpose_mutation_loci(mutation_loci, cssplits_sample)
        cssplits_mutation_loci_sample = _extract_cssplits_in_mutation_by_3mer(cssplits_sample, mutation_loci_transposed)
        cssplits_mutation_loci_control = _extract_cssplits_in_mutation_by_3mer(cssplits_control, mutation_loci_transposed)
        mutation_score = make_score(cssplits_mutation_loci_sample, cssplits_mutation_loci_control, knockin_loci)
        scores_sample = _annotate_score(cssplits_mutation_loci_sample, mutation_score)
        scores_control = _annotate_score(cssplits_mutation_loci_control, mutation_score)
        labels = return_labels(scores_sample, scores_control)
        labels_reorder = _reorder_labels(labels, start=max_label)
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
