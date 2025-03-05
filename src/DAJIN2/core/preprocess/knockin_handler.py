from __future__ import annotations

import pickle
from pathlib import Path

import midsv

from DAJIN2.core.preprocess import mapping

###########################################################
# Consider all mutations are possible in the knockin region
# For large deletion alleles, the deleted sequence becomes the knock-in region, so all mutations within this region are taken into consideration.
# The code identifies the flox knock-in sites as deletions not present in the control.
###########################################################


def is_valid_file(control_fasta: Path, knockin_fasta: Path) -> bool:
    return knockin_fasta != control_fasta and knockin_fasta.suffix == ".fasta"


def get_index_of_knockin_loci(control_fasta: Path, knockin_fasta: Path) -> set[int]:
    alignments = mapping.to_sam(knockin_fasta, control_fasta, preset="map-ont")
    alignments = [a.split("\t") for a in alignments]
    alignments_midsv = list(midsv.transform(alignments, midsv=False, cssplit=True, qscore=False))
    if not alignments_midsv:
        return set()
    alignments_midsv = next(iter(midsv.transform(alignments, midsv=False, cssplit=True, qscore=False)))
    cssplits = alignments_midsv["CSSPLIT"].split(",")

    return {i for i, cs in enumerate(cssplits) if cs.startswith("-")}


def extract_knockin_loci(TEMPDIR: str | Path, SAMPLE_NAME: str) -> None:
    control_fasta = Path(TEMPDIR, SAMPLE_NAME, "fasta", "control.fasta")
    for knockin_fasta in Path(TEMPDIR, SAMPLE_NAME, "fasta").iterdir():
        if not is_valid_file(control_fasta, knockin_fasta):
            continue

        allele = knockin_fasta.stem
        path_output = Path(TEMPDIR, SAMPLE_NAME, "knockin_loci", allele, "knockin.pickle")
        path_output.parent.mkdir(parents=True, exist_ok=True)
        if path_output.exists():
            continue

        knockin_loci = get_index_of_knockin_loci(control_fasta, knockin_fasta)

        with open(path_output, "wb") as p:
            pickle.dump(knockin_loci, p)
