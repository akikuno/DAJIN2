from __future__ import annotations

import sys
from itertools import groupby
from pathlib import Path

import pandas as pd
import wslPath

from DAJIN2.core import core_execute
from DAJIN2.postprocess import report
from DAJIN2.preprocess.validate_inputs import (validate_files,
                                               validate_genome_and_fetch_urls)


def execute_batch_mode(arguments: dict[str]):
    threads = int(arguments["threads"])
    path_batchfile = arguments["file"]

    ###############################################################################
    # Validate batch file
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

    # ----------------------------------------------------------
    # Load the file
    # ----------------------------------------------------------
    try:
        df_batchfile = pd.read_excel(path_batchfile)
        inputs = []
        inputs.append(df_batchfile.columns.to_list())
        inputs += df_batchfile.values.tolist()
    except ValueError:
        inputs = [s.split(",") for s in Path(path_batchfile).read_text().strip().split("\n")]

    ################################################################################
    # Validate Column of the batch file
    ################################################################################
    columns = inputs[0]
    columns_set = set(columns)
    columns_required = set(["sample", "control", "allele", "name"])
    columns_all = set(["sample", "control", "allele", "name", "genome"])

    if not columns_required.issubset(columns_set):
        raise ValueError(f"{path_batchfile} must contain 'sample', 'control', and 'allele' in the header")

    if not columns_set.issubset(columns_all):
        raise ValueError(
            f"Accepted header names of {path_batchfile} are 'sample', 'control', 'allele', 'name', and 'genome'."
        )

    ################################################################################
    # Validate contents in the batch file
    ################################################################################
    index_name = columns.index("name")
    contents = inputs[1:]
    contents.sort(key=lambda x: x[index_name])

    for name, groups in groupby(contents, key=lambda x: x[index_name]):
        for group in groups:
            args = {h: g for h, g in zip(columns, group)}
            validate_files(args["sample"], args["control"], args["allele"])
            if "genome" in args:
                UCSC_URL, GOLDENPATH_URL = validate_genome_and_fetch_urls(args["genome"])
            else:
                UCSC_URL, GOLDENPATH_URL = None, None
    ##############################################################################
    # Perform DAJIN
    ##############################################################################

    for name, groups in groupby(contents, key=lambda x: x[index_name]):
        done_controle = False
        for group in groups:
            args = {h: g for h, g in zip(columns, group)}
            args["threads"] = threads
            if done_controle == False:
                print(f"{args['control']} is now processing...", file=sys.stderr)
                core_execute.execute_control(args)
                done_controle = True
            print(f"{args['sample']} is now processing...", file=sys.stderr)
            core_execute.execute_sample(args)
        report.report(name)
        print(f"Finished! Open DAJINResults/{name} to see the report.", file=sys.stderr)
