import midsv
from pathlib import Path


def output_midsv(TEMPDIR, path_sam, DICT_ALLELE):
    name_fasta = path_sam.stem.split("_")[1]
    if name_fasta not in set(DICT_ALLELE.keys()):
        return
    sam = midsv.read_sam(path_sam)
    midsv_jsonl = midsv.transform(sam, midsv=False, cssplit=True, qscore=False)
    output_jsonl = Path(TEMPDIR, "midsv", f"{path_sam.stem}.jsonl")
    midsv.write_jsonl(midsv_jsonl, output_jsonl)
