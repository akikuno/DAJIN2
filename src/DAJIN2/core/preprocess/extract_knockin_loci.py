from __future__ import annotations

from pathlib import Path
import pickle
import midsv

from DAJIN2.core.preprocess import mappy_align


###########################################################
# Consider all mutations are possible in the knockin region
###########################################################


def extract_knockin_loci(TEMPDIR: str | Path, SAMPLE_NAME: str) -> None:
    reference = Path(TEMPDIR, "fasta", "control.fasta")
    for query in Path(TEMPDIR, "fasta").iterdir():
        if query.suffix != ".fasta":
            continue
        if reference == query:
            continue
        allele = query.stem
        if Path(TEMPDIR, "knockin_loci", f"{SAMPLE_NAME}_{allele}.pickle").exists():
            continue
        alignments = mappy_align.to_sam(reference, query, preset="splice")
        alignments = [a.split("\t") for a in alignments]
        alignments_midsv = midsv.transform(alignments, midsv=False, cssplit=True, qscore=False)[0]
        cssplits = alignments_midsv["CSSPLIT"].split(",")
        knockin_loci = {i for i, cs in enumerate(cssplits) if cs == "N" or cs.startswith("-")}
        with open(Path(TEMPDIR, "knockin_loci", f"{SAMPLE_NAME}_{allele}.pickle"), "wb") as p:
            pickle.dump(knockin_loci, p)
