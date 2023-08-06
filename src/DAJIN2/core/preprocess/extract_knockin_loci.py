from __future__ import annotations

from pathlib import Path
import pickle
import midsv

from DAJIN2.core.preprocess import align


###########################################################
# Consider all mutations are possible in the knockin region
###########################################################


def extract_knockin_loci(TEMPDIR: str | Path, SAMPLE_NAME: str) -> None:
    reference = Path(TEMPDIR, SAMPLE_NAME, "fasta", "control.fasta")
    for query in Path(TEMPDIR, SAMPLE_NAME, "fasta").iterdir():
        if query.suffix != ".fasta":
            continue
        if reference == query:
            continue
        allele = query.stem
        path_output = Path(TEMPDIR, SAMPLE_NAME, "knockin_loci", f"{allele}.pickle")
        if path_output.exists():
            continue
        alignments = align.to_sam(reference, query, preset="splice")
        alignments = [a.split("\t") for a in alignments]
        alignments_midsv = midsv.transform(alignments, midsv=False, cssplit=True, qscore=False)[0]
        cssplits = alignments_midsv["CSSPLIT"].split(",")
        knockin_loci = {i for i, cs in enumerate(cssplits) if cs == "N" or cs.startswith("-")}
        with open(path_output, "wb") as p:
            pickle.dump(knockin_loci, p)
