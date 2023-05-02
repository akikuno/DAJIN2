from __future__ import annotations
from collections import defaultdict
from pathlib import Path

import midsv


def _calc_percent_indels(cssplits):
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


def _extract(cssplits_sample, cssplits_control) -> dict[int, str]:
    percent_sample = _calc_percent_indels(cssplits_sample)
    percent_control = _calc_percent_indels(cssplits_control)
    mutation_loci = dict()
    for i, (samp, cont) in enumerate(zip(percent_sample, percent_control)):
        for mutation_type in ["ins", "del", "sub"]:
            if abs(samp[mutation_type] - cont[mutation_type]) > 0.5:
                mutation_loci.update({i: mutation_type})
    return mutation_loci


###############################################################################


def extract_mutation_loci(TEMPDIR: Path, FASTA_ALLELES: dict, SAMPLE_NAME: str, CONTROL_NAME: str) -> dict[str, dict[int, str]]:
    mutation_loci = dict()
    for allele in FASTA_ALLELES.keys():
        midsv_sample = midsv.read_jsonl((Path(TEMPDIR, "midsv", f"{SAMPLE_NAME}_{allele}.jsonl")))
        midsv_control = midsv.read_jsonl((Path(TEMPDIR, "midsv", f"{CONTROL_NAME}_{allele}.jsonl")))
        cssplits_sample = [cs["CSSPLIT"].split(",") for cs in midsv_sample]
        cssplits_control = [cs["CSSPLIT"].split(",") for cs in midsv_control]
        mutation_loci.update({allele: _extract(cssplits_sample, cssplits_control)})
    return mutation_loci
