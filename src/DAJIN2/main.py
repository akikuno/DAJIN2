from __future__ import annotations

import os
import argparse
import multiprocessing
import sys
from itertools import groupby, islice
from pathlib import Path

import pandas as pd
import wslPath

from DAJIN2 import gui, view
from DAJIN2.core import core_execute
from DAJIN2.postprocess import report
from DAJIN2.preprocess.validate_inputs import validate_files, validate_genome_and_fetch_urls

VERSION = "0.2.3"

# prevent BLAS from using all cores
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"


def execute_single_mode(arguments: dict[str]):
    ################################################################################
    # Validate contents
    ################################################################################
    validate_files(arguments["sample"], arguments["control"], arguments["allele"])
    URLS_GENOME: dict[str, dict[str, str]] = dict()
    if "genome" in arguments:
        URLS_GENOME.update({arguments["genome"]: validate_genome_and_fetch_urls(arguments["genome"])})
        arguments.update(URLS_GENOME[arguments["genome"]])
    ##############################################################################
    # Perform DAJIN2
    ##############################################################################
    core_execute.execute_control(arguments)
    core_execute.execute_sample(arguments)
    name = arguments["name"]
    report.report(name)
    print(
        f"\N{party popper} Finished! Open DAJINResults/{name} to see the report.",
        file=sys.stderr,
    )


def _update_threads(threads) -> int:
    threads_updated = min(int(threads), os.cpu_count() - 1)
    threads_updated = max(1, threads_updated)
    return threads_updated


def _extract_unique_contents(list: list) -> list:
    list_unique = []
    for content in list:
        if content not in list_unique:
            list_unique.append(content)
    return list_unique


def _batched(iterable, chunk_size):
    iterator = iter(iterable)
    chunk = tuple(islice(iterator, chunk_size))
    while chunk:
        yield chunk
        chunk = tuple(islice(iterator, chunk_size))


def _run_multiprocess(function, arguments: list, num_workers: int = 1) -> None:
    arguments_batched = _batched(arguments, num_workers)
    for args in arguments_batched:
        jobs = []
        for arg in args:
            p = multiprocessing.Process(target=function, args=(arg,))
            jobs.append(p)
            p.start()
        for job in jobs:
            if job.exitcode == 1:
                sys.exit(1)
            job.join()
    return


def execute_batch_mode(arguments: dict[str]):
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
    if not {"sample", "control", "allele", "name"}.issubset(set(columns)):
        raise ValueError(f"{path_batchfile} must contain 'sample', 'control', and 'allele' in the header")
    if not set(columns).issubset({"sample", "control", "allele", "name", "genome"}):
        raise ValueError(
            f"Accepted header names of {path_batchfile} are 'sample', 'control', 'allele', 'name', and 'genome'."
        )
    ################################################################################
    # Validate contents in the batch file
    ################################################################################
    index_name = columns.index("name")
    contents = inputs[1:]
    contents.sort(key=lambda x: x[index_name])
    for _, groups in groupby(contents, key=lambda x: x[index_name]):
        for group in groups:
            args = {h: g for h, g in zip(columns, group)}
            validate_files(args["sample"], args["control"], args["allele"])
            URLS_GENOME: dict[str, dict[str, str]] = dict()
            if "genome" in args and not args["genome"] in URLS_GENOME:
                URLS_GENOME.update({args["genome"]: validate_genome_and_fetch_urls(args["genome"])})
    ##############################################################################
    # Perform DAJIN2
    ##############################################################################
    num_workers = _update_threads(arguments["threads"])
    index_name = columns.index("name")
    contents.sort(key=lambda x: x[index_name])
    contents_sample = []
    for name, groups in groupby(contents, key=lambda x: x[index_name]):
        groups = list(groups)
        # ------------------------------
        # Handle controls
        # ------------------------------
        contents_control = []
        for group in groups:
            args = {h: g for h, g in zip(columns, group)}
            args["sample"] = args["control"]
            if len(groups) > 1:
                args["threads"] = 1
            if "genome" in args:
                args.update(URLS_GENOME[args["genome"]])
            contents_control.append(args)
        contents_control_unique = _extract_unique_contents(contents_control)
        _run_multiprocess(core_execute.execute_control, contents_control_unique, num_workers)
        # ------------------------------
        # Handle samples
        # ------------------------------
        contents_sample = []
        for group in groups:
            args = {h: g for h, g in zip(columns, group)}
            if args["sample"] == args["control"]:
                continue
            if "threads" not in args:
                args["threads"] = 1
            args.update(URLS_GENOME[args["genome"]])
            contents_sample.append(args)
        contents_sample_unique = _extract_unique_contents(contents_sample)
        _run_multiprocess(core_execute.execute_sample, contents_sample_unique, num_workers)
        report.report(name)
        print(
            f"\N{party popper} Finished! Open DAJINResults/{name} to see the report.",
            file=sys.stderr,
        )


def execute():
    parser = argparse.ArgumentParser()

    ###############################################################################
    # Single mode
    ###############################################################################

    parser.add_argument("-s", "--sample", type=str, help="Full path to a sample FASTQ file")
    parser.add_argument("-c", "--control", type=str, help="Full path to a control FASTQ file")
    parser.add_argument("-a", "--allele", type=str, help="Full path to a FASTA file")
    parser.add_argument("-n", "--name", type=str, help="Output directory name")
    parser.add_argument(
        "-g", "--genome", type=str, default="", help="Reference genome ID (e.g hg38, mm10) [default: '']"
    )
    parser.add_argument("-t", "--threads", type=int, default=1, help="Number of threads [default: 1]")
    parser.add_argument("-v", "--version", action="version", version=f"DAJIN2 version {VERSION}")
    parser.add_argument("-d", "--debug", action="store_true", help=argparse.SUPPRESS)

    ###############################################################################
    # Batch mode
    ###############################################################################

    def batchmode(args):
        arguments = dict()
        arguments["file"] = args.file
        arguments["threads"] = int(args.threads)
        arguments["debug"] = args.debug
        execute_batch_mode(arguments)

    subparser = parser.add_subparsers()
    parser_batch = subparser.add_parser("batch", help="DAIJN2 batch mode")
    parser_batch.add_argument("-f", "--file", required=True, type=str, help="CSV or Excel file.")
    parser_batch.add_argument("-t", "--threads", default=1, type=int, help="Number of threads [default: 1]")
    parser_batch.add_argument("-d", "--debug", action="store_true", help=argparse.SUPPRESS)
    parser_batch.set_defaults(handler=batchmode)

    ###############################################################################
    # GUI mode
    ###############################################################################

    def guimode(args):
        gui.execute()

    parser_gui = subparser.add_parser("gui", help="DAIJN2 GUI mode")
    parser_gui.set_defaults(handler=guimode)

    ###############################################################################
    # View mode
    ###############################################################################

    def viewmode(args):
        view.execute(args.name)

    parser_view = subparser.add_parser("view", help="DAIJN2 View mode to launch igvjs")
    parser_view.add_argument("-n", "--name", required=True, type=str, help="Output name of the report")
    parser_view.set_defaults(handler=viewmode)

    ###############################################################################
    # Parse arguments
    ###############################################################################

    args = parser.parse_args()

    if hasattr(args, "handler"):
        args.handler(args)
    else:
        if args.sample is None:
            raise AttributeError("the following arguments are required: -s/--sample")
        if args.control is None:
            raise AttributeError("the following arguments are required: -c/--control")
        if args.allele is None:
            raise AttributeError("the following arguments are required: -a/--allele")
        if args.name is None:
            raise AttributeError("the following arguments are required: -n/--name")
        arguments = dict()
        arguments["sample"] = args.sample
        arguments["control"] = args.control
        arguments["allele"] = args.allele
        arguments["name"] = args.name
        if args.genome:
            arguments["genome"] = args.genome
        arguments["threads"] = int(args.threads)
        arguments["debug"] = args.debug
        execute_single_mode(arguments)


if __name__ == "__main__":
    execute()
