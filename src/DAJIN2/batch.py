from __future__ import annotations
from pathlib import Path
from itertools import groupby
import pandas as pd
import wslPath
from .core import core
from .postprocess import report


def execute(arguments: dict[str]):
    path_batchfile = arguments["file"]
    threads = arguments["threads"]
    # debug = arguments["debug"]

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
        inputs.append(df_batchfile.columns.to_list())
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
        for group in groups:
            args = {h: g for h, g in zip(keys, group)}
            args["threads"] = threads
            core.execute(args)
        report.report(name)
