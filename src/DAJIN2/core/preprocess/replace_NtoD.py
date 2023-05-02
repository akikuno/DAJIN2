from __future__ import annotations

from pathlib import Path

import midsv


def _replaceNtoD(cssplits_sample, sequence) -> list[list[str]]:
    cssplits_replaced = cssplits_sample.copy()
    for i, cssplits in enumerate(cssplits_sample):
        # extract right/left index of the end of sequential Ns
        left_idx_n = 0
        for cs in cssplits:
            if cs != "N":
                break
            left_idx_n += 1
        right_idx_n = 0
        for cs in cssplits[::-1]:
            if cs != "N":
                break
            right_idx_n += 1
        right_idx_n = len(cssplits) - right_idx_n - 1
        # replace sequential Ns within the sequence
        for j, (cs, seq) in enumerate(zip(cssplits, sequence)):
            if left_idx_n <= j <= right_idx_n and cs == "N":
                cssplits_replaced[i][j] = f"-{seq}"
    return cssplits_replaced


###############################################################################
# main
###############################################################################

def replace_NtoD(midsv_sample_alleles: dict[str, list[dict]], FASTA_ALLELES: dict) -> dict[str, list[dict]]:
    """
    Convert any `N` as deletions other than consecutive `N` from both ends
    """
    midsv_updated = dict()
    for allele, sequence in FASTA_ALLELES.items():
        midsv_sample = midsv_sample_alleles[allele]
        cssplits_sample = [cs["CSSPLIT"].split(",") for cs in midsv_sample]
        cssplits_replaced = _replaceNtoD(cssplits_sample, sequence)
        midsv_cssplits = [",".join(cs) for cs in cssplits_replaced]
        # Save as a json
        for i, cssplits in enumerate(midsv_cssplits):
            midsv_sample[i]["CSSPLIT"] = cssplits
        midsv_updated[allele] = midsv_sample
    return midsv_updated
