from pathlib import Path
import midsv
from collections import defaultdict


def mask_control(TEMPDIR: Path, DICT_ALLELE: dict, CONTROL_NAME: str):
    for allele in DICT_ALLELE.keys():
        if allele == "control":
            continue
        masks_control = defaultdict(list)
        path_midsv = Path(TEMPDIR, "midsv", f"{CONTROL_NAME}_{allele}.jsonl")
        cssplit_control = midsv.read_jsonl(path_midsv)
        cssplit_control = [cs["CSSPLIT"].split(",") for cs in cssplit_control]
        coverage = len(cssplit_control)
        cssplit_control_transposed = list(zip(*cssplit_control))
        masks_loci = [False] * (len(cssplit_control_transposed) - 1)
        for i, cs in enumerate(cssplit_control_transposed):
            count_match = sum(1 for c in cs if c.startswith("=") or c == "N")
            percent_match = count_match / coverage * 100
            if percent_match < 1:
                masks_loci[i] = True
        masks_control[allele] = masks_loci
    return masks_control
