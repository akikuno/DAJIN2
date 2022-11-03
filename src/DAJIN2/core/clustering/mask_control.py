from collections import defaultdict
from pathlib import Path

import midsv


def mask_control(TEMPDIR: Path, DICT_ALLELE: dict, CONTROL_NAME: str):
    """
    Mask loci that have less than 1% of matches in control seqence.
    Alignments between control and knock-in alleles produses deletion loci in control.
    The deletion loci are not sequencing error, so they should be ignored.
    """
    masks_control = defaultdict(list)
    for allele in DICT_ALLELE.keys():
        if allele == "control":
            continue
        path_midsv = Path(TEMPDIR, "midsv", f"{CONTROL_NAME}_{allele}.jsonl")
        cssplit_control = midsv.read_jsonl(path_midsv)
        cssplit_control = [cs["CSSPLIT"].split(",") for cs in cssplit_control]
        coverage = len(cssplit_control)
        cssplit_control_transposed = list(zip(*cssplit_control))
        masks_loci = []
        for i, cs in enumerate(cssplit_control_transposed):
            count_match = sum(1 for c in cs if c.startswith("="))
            percent_match = count_match / coverage * 100
            if percent_match < 1:
                masks_loci.append(True)
            else:
                masks_loci.append(False)
        masks_control[allele] = masks_loci
    return masks_control
