from __future__ import annotations

from collections import defaultdict
from itertools import groupby
from typing import Generator
import json
from pathlib import Path
from DAJIN2.core.clustering.make_score import make_score
from DAJIN2.core.clustering.return_labels import return_labels
import random

# def _compress_insertion(cssplits: Generator[list[str]]) -> Generator[dict[str, int]]:
#     """Insertion will be subdivided by sequence error in the its sequence,
#     so it is compressed as a '+I' to eliminate mutations.
#     #TODO ただ、これでは、insertion配列の中に真のmutationがある場合に
#     #TODO そのmutationを抽出できないので**insertion配列の中にmutationがある場合は
#     #TODO insertion配列をそのまま残す**必要がある。
#     """
#     for cssplit in cssplits:
#         for i, cs in enumerate(cssplit):
#             if cs.startswith("+"):
#                 cssplit[i] = "+I" + cs.split("|")[-1]
#         yield cssplit


def _generate_mutation_kmers(midsv_sample: Generator[list[str]], mutation_loci: list[set[str]]) -> Generator[list[str]]:
    for cssplit in (cs["CSSPLIT"].split(",") for cs in midsv_sample):
        cs_mutation = ["N,N,N"]
        for i in range(1, len(cssplit) - 1):
            if mutation_loci[i] == set():
                cs_mutation.append("N,N,N")
                continue
            mutation = cssplit[i][0]  # +, - , *, =, N
            """Insertion will be subdivided by sequence error in the its sequence,
            so it is compressed as a '+I' to eliminate mutations.
            #TODO ただ、これでは、insertion配列の中に真のmutationがある場合に
            #TODO そのmutationを抽出できないので**insertion配列の中にmutationがある場合は
            #TODO insertion配列をそのまま残す**必要がある。
            """
            if mutation == "+":
                cssplit[i] = "+I" + cssplit[i].split("|")[-1]
            if mutation in mutation_loci[i]:
                kmer = ",".join([cssplit[i - 1], cssplit[i], cssplit[i + 1]])
                cs_mutation.append(kmer)
            else:
                cs_mutation.append("N,N,N")
        cs_mutation.append("N,N,N")
        yield cs_mutation


def _annotate_score(cssplits: Generator[list[str]], mutation_score: list[dict[str:float]]) -> Generator[list[float]]:
    for cssplit in cssplits:
        score = [0 for _ in range(len(cssplit))]
        for i, (cs, mutscore) in enumerate(zip(cssplit, mutation_score)):
            if mutscore == {}:
                continue
            if cs in mutscore:
                score[i] = mutscore[cs]
        yield score


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


def read_json(filepath: Path | str) -> Generator[dict[str, str]]:
    with open(filepath, "r") as f:
        for line in f:
            yield json.loads(line)


def write_json(filepath: Path | str, data: Generator) -> None:
    with open(filepath, "w") as f:
        for line in data:
            f.write(json.dumps(line) + "\n")


def add_labels(
    classif_sample, TEMPDIR, SAMPLE_NAME, CONTROL_NAME, MUTATION_LOCI_ALLELES, KNOCKIN_LOCI_ALLELES, THREADS: int = 1
) -> list[dict[str]]:
    labels_all = []
    max_label = 0
    classif_sample.sort(key=lambda x: x["ALLELE"])
    for allele, group in groupby(classif_sample, key=lambda x: x["ALLELE"]):
        RANDOM_NUM = random.randint(0, 10**10)
        mutation_loci: dict[str, set[str]] = MUTATION_LOCI_ALLELES[allele]
        if all(m == set() for m in mutation_loci):
            labels = [1] * len(classif_sample)
            labels_reorder = _reorder_labels(labels, start=max_label)
            max_label = max(labels_reorder)
            labels_all.extend(labels_reorder)
            continue
        if allele in KNOCKIN_LOCI_ALLELES:
            knockin_loci = KNOCKIN_LOCI_ALLELES[allele]
        else:
            knockin_loci = set()
        path_sample = Path(TEMPDIR, "clustering", f"{SAMPLE_NAME}_{allele}_{RANDOM_NUM}.json")
        path_control = Path(TEMPDIR, "midsv", f"{CONTROL_NAME}_{allele}.json")
        write_json(path_sample, group)
        mutation_score = make_score(
            _generate_mutation_kmers(read_json(path_sample), mutation_loci),
            _generate_mutation_kmers(read_json(path_control), mutation_loci),
            mutation_loci,
            knockin_loci,
        )
        scores_sample = _annotate_score(_generate_mutation_kmers(read_json(path_sample), mutation_loci), mutation_score)
        scores_control = _annotate_score(
            _generate_mutation_kmers(read_json(path_control), mutation_loci), mutation_score
        )
        path_score_sample = Path(TEMPDIR, "clustering", f"{SAMPLE_NAME}_{allele}_score_{RANDOM_NUM}.json")
        path_score_control = Path(TEMPDIR, "clustering", f"{CONTROL_NAME}_{allele}_score_{RANDOM_NUM}.json")
        write_json(path_score_sample, scores_sample)
        write_json(path_score_control, scores_control)
        labels = return_labels(path_score_sample, path_score_control)
        labels_reorder = _reorder_labels(labels, start=max_label)
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
