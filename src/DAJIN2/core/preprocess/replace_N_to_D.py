from __future__ import annotations
from copy import deepcopy
from pathlib import Path

import midsv


def replaceNtoD(cssplits_sample, sequence):
    cssplits_replaced = deepcopy(cssplits_sample)
    for i, cssplits in enumerate(cssplits_sample):
        flag_n_start = True
        flag_n_end = True
        for j, (start, end) in enumerate(zip(cssplits, cssplits[::-1])):
            if j == (len(cssplits) + 1) // 2:
                break
            if flag_n_start and start != "N":
                flag_n_start = False
            if flag_n_end and end != "N":
                flag_n_end = False
            if not flag_n_start and start == "N":
                cssplits_replaced[i][j] = f"-{sequence[j]}"
            if not flag_n_end and end == "N":
                j_inv = len(cssplits) - j - 1
                cssplits_replaced[i][j_inv] = f"-{sequence[j_inv]}"
    return cssplits_replaced


###############################################################################
# main
###############################################################################


def execute(TEMPDIR: Path, FASTA_ALLELES: dict, SAMPLE_NAME: str) -> None:
    """
    Convert any `N` as deletions other than consecutive `N` from both ends
    """
    for allele, sequence in FASTA_ALLELES.items():
        midsv_sample = midsv.read_jsonl((Path(TEMPDIR, "midsv", f"{SAMPLE_NAME}_splice_{allele}.jsonl")))
        cssplits_sample = [cs["CSSPLIT"].split(",") for cs in midsv_sample]
        cssplits_replaced = replaceNtoD(cssplits_sample, sequence)
        midsv_cssplits = [",".join(cs) for cs in cssplits_replaced]
        # Save as a json
        for i, cssplits in enumerate(midsv_cssplits):
            midsv_sample[i]["CSSPLIT"] = cssplits
        midsv.write_jsonl(midsv_sample, Path(TEMPDIR, "midsv", f"{SAMPLE_NAME}_splice_{allele}.jsonl"))
