from pathlib import Path

import midsv


def output_midsv(TEMPDIR, path_sam):
    sam = midsv.read_sam(path_sam)
    midsv_jsonl = midsv.transform(sam, midsv=False, cssplit=True, qscore=True)
    output_jsonl = Path(TEMPDIR, "midsv", f"{path_sam.stem}.jsonl")
    midsv.write_jsonl(midsv_jsonl, output_jsonl)
