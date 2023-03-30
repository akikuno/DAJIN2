from __future__ import annotations
import midsv
from pathlib import Path
from collections import defaultdict


def calc_percent_indels(cssplits):
    percelt_indels = []
    cssplits_transposed = [list(t) for t in zip(*cssplits)]
    for cssplit in cssplits_transposed:
        coverage = 0
        count_indelsub = {"ins": 0, "del": 0, "sub": 0}
        for cs in cssplit:
            if cs == "N":
                continue
            coverage += 1
            if cs.startswith("+"):
                count_indelsub["ins"] += 1
            elif cs.startswith("-"):
                count_indelsub["del"] += 1
            elif cs.startswith("*"):
                count_indelsub["sub"] += 1
        if coverage == 0:
            per_indels = {"ins": 0, "del": 0, "sub": 0}
        else:
            per_indels = {mutation: (count / coverage * 100) for mutation, count in count_indelsub.items()}
        percelt_indels.append(per_indels)
    return percelt_indels


def extract(cssplits_sample, cssplits_control) -> set():
    percent_sample = calc_percent_indels(cssplits_sample)
    percent_control = calc_percent_indels(cssplits_control)
    mutation_loci = set()
    for i, (samp, cont) in enumerate(zip(percent_sample, percent_control)):
        for mutation_type in ["ins", "del", "sub"]:
            if i in mutation_loci:
                break
            if abs(samp[mutation_type] - cont[mutation_type]) > 1:
                mutation_loci.add(i)
    return mutation_loci


###############################################################################


def main(TEMPDIR: Path, FASTA_ALLELES: dict, SAMPLE_NAME: str, CONTROL_NAME:str) -> defaultdict(set):
    mutation_loci = defaultdict(set)
    for allele in FASTA_ALLELES.keys():
        midsv_sample = midsv.read_jsonl((Path(TEMPDIR, "midsv", f"{SAMPLE_NAME}_splice_{allele}.jsonl")))
        midsv_control = midsv.read_jsonl((Path(TEMPDIR, "midsv", f"{CONTROL_NAME}_splice_{allele}.jsonl")))
        cssplits_sample = [cs["CSSPLIT"].split(",") for cs in midsv_sample]
        cssplits_control = [cs["CSSPLIT"].split(",") for cs in midsv_control]
        mutation_loci[allele] = extract(cssplits_sample, cssplits_control)
    return mutation_loci
