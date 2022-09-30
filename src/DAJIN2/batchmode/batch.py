from __future__ import annotations
from pathlib import Path
from itertools import groupby
from ..core import main
import wslPath
import midsv
import pandas as pd


def batch(arguments: dict[str]):
    path_batchfile = arguments["file"]
    threads = arguments["threads"]
    debug = arguments["debug"]

    ###############################################################################
    # Check arguments
    ###############################################################################
    try:
        path_batchfile = wslPath.toPosix(path_batchfile)
    except ValueError:
        pass

    # ----------------------------------------------------------
    # Check file exists
    # ----------------------------------------------------------
    if not Path(path_batchfile).exists():
        raise FileNotFoundError(f"'{path_batchfile}' does not exist.")

    try:
        df_batchfile = pd.read_excel(path_batchfile)
        inputs = []
        inputs = df_batchfile.columns.to_list()
        inputs += df_batchfile.values.tolist()
    except ValueError:
        inputs = [s.split(",") for s in Path(path_batchfile).read_text().strip().split("\n")]

    # ----------------------------------------------------------
    # Check Column
    # ----------------------------------------------------------
    keys = inputs[0]
    keys_set = set(keys)
    keys_required = set(["sample", "control", "allele", "name"])
    keys_all = set(["sample", "control", "allele", "name", "genome"])

    if not keys_required.issubset(keys_set):
        raise ValueError(f"{path_batchfile} must contain 'sample', 'control', and 'allele' in the header")

    if not keys_set.issubset(keys_all):
        raise ValueError(
            f"Accepted header names of {path_batchfile} are 'sample', 'control', 'allele', 'name', and 'genome'."
        )

    ##############################################################################
    # Perform DAJIN
    ##############################################################################
    index_name = [i for i, key in enumerate(keys) if key == "name"][0]
    contents = inputs[1:]
    contents.sort(key=lambda x: x[index_name])

    for name, groups in groupby(contents, key=lambda x: x[index_name]):
        output_dir = Path("DAJINResults", f".tempdir_{name}", "batch")
        output_dir.mkdir(exist_ok=True, parents=True)
        for group in groups:
            args = {h: g for h, g in zip(keys, group)}
            args["threads"] = threads
            SAMPLE_NAME, RESULT = main.main(args)
            # output the result
            midsv.write_jsonl(RESULT, Path(output_dir, f"{SAMPLE_NAME}.jsonl"))


##############################################################################
# Reports
##############################################################################
name = "Ayabe-Task1"
output_dir = Path("DAJINResults", f".tempdir_{name}", "batch")
df_results = pd.DataFrame()
for path_result in output_dir.iterdir():
    SAMPLE_NAME = path_result.stem
    result = midsv.read_jsonl(path_result)
    df_result = pd.DataFrame(result)
    df_result["SAMPLE"] = SAMPLE_NAME
    df_results = pd.concat([df_results, df_result])

