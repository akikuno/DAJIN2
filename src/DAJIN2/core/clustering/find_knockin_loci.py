from __future__ import annotations
from collections import defaultdict, Counter
from pathlib import Path
from copy import deepcopy
import midsv


def replaceNtoMatch(cssplit: list[str], sequence) -> list[str]:
    cssplit_replaced = deepcopy(cssplit)
    # Replace N to @ at the left ends
    for i, cs in enumerate(cssplit_replaced):
        if cs != "N":
            break
        cssplit_replaced[i] = "=" + sequence[i]
    # Replace N to @ at the right ends
    cssplit_replaced = cssplit_replaced[::-1]
    for i, cs in enumerate(cssplit_replaced):
        if cs != "N":
            break
        cssplit_replaced[i] = "=" + sequence[::-1][i]
    cssplit_replaced = cssplit_replaced[::-1]
    return cssplit_replaced


def call_consensus_80percent(cssplit, sequence) -> list[str]:
    coverage = len(cssplit)
    cssplit_transposed = list(zip(*cssplit))
    cssplit_mostcommon = [Counter(cs).most_common()[0] for cs in cssplit_transposed]
    cssplit_consensus = []
    for i, (key, val) in enumerate(cssplit_mostcommon):
        val_percent = val / coverage * 100
        if key.startswith("="):
            cssplit_consensus.append(key)
        elif val_percent > 80:
            cssplit_consensus.append(key)
        else:
            # Mutations occuring in less than 80% are considered sequencing errors
            cssplit_consensus.append(sequence[i])
    return cssplit_consensus


def find_knockin_loci(TEMPDIR: Path, DICT_ALLELE: dict, CONTROL_NAME: str) -> defaultdict(set[int]):
    """
    Alignments between control and knock-in alleles produses deletion loci in control.
    The deletion loci are knock-in loci, not sequencing error, so they should be ignored.
    """
    knockin_loci = defaultdict(set[int])
    for allele in DICT_ALLELE.keys():
        if allele == "control":
            continue
        path_midsv = Path(TEMPDIR, "midsv", f"{CONTROL_NAME}_{allele}.jsonl")
        sequence = DICT_ALLELE[allele]
        cssplit = midsv.read_jsonl(path_midsv)
        cssplit = [cs["CSSPLIT"].split(",") for cs in cssplit]
        cssplit = [replaceNtoMatch(cs, sequence) for cs in cssplit]
        cssplit_consensus = call_consensus_80percent(cssplit, sequence)
        knockin_loci[allele] = {i for i, cs in enumerate(cssplit_consensus) if cs.startswith("-")}
    return knockin_loci
