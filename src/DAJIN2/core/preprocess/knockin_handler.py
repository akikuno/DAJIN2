from __future__ import annotations

import pickle
from pathlib import Path

import midsv

from DAJIN2.core.preprocess import mapping


###########################################################
# Consider all mutations are possible in the knockin region
# For large deletion alleles, the deleted sequence becomes the knock-in region, so all mutations within this region are taken into consideration.
# TODO: Maybe this knock-in consideration isn't necessary?
###########################################################


def _is_invalid_file(reference: Path, query: Path) -> bool:
    return query == reference or query.suffix != ".fasta"


def _get_knockin_loci(reference: Path, query: Path) -> set[int]:
    alignments = mapping.to_sam(reference, query, preset="splice")
    alignments = [a.split("\t") for a in alignments]
    alignments_midsv = midsv.transform(alignments, midsv=False, cssplit=True, qscore=False)[0]
    cssplits = alignments_midsv["CSSPLIT"].split(",")

    return {i for i, cs in enumerate(cssplits) if cs.startswith("-")}


def extract_knockin_loci(TEMPDIR: str | Path, SAMPLE_NAME: str) -> None:
    reference = Path(TEMPDIR, SAMPLE_NAME, "fasta", "control.fasta")
    for query in Path(TEMPDIR, SAMPLE_NAME, "fasta").iterdir():
        if _is_invalid_file(reference, query):
            continue

        allele = query.stem
        path_output = Path(TEMPDIR, SAMPLE_NAME, "knockin_loci", f"{allele}.pickle")
        if path_output.exists():
            continue

        knockin_loci = _get_knockin_loci(reference, query)

        with open(path_output, "wb") as p:
            pickle.dump(knockin_loci, p)
