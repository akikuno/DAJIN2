from __future__ import annotations

import os
import sys
import logging
import datetime
import argparse

import traceback
from multiprocessing import Process, Queue

from pathlib import Path
from itertools import groupby, islice
from typing import Generator

from DAJIN2 import gui, view
from DAJIN2.core import core
from DAJIN2.utils import config, io, report_generator, input_validator


VERSION = "0.3.1"

# prevent BLAS from using all cores
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"

###########################################################
# Setting logger
###########################################################

current_time = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")

stderr_handler = logging.StreamHandler()
stderr_handler.setFormatter(logging.Formatter("%(message)s"))

file_handler = logging.FileHandler(f"{current_time}_DAJIN2.log")
file_handler.setFormatter(logging.Formatter("%(message)s"))

logging.basicConfig(
    level=logging.INFO,
    format="[%(asctime)s] [%(levelname)s] [%(processName)s] %(message)s",
    handlers=[stderr_handler, file_handler],
)


def update_threads(threads: int) -> int:
    available_threads = os.cpu_count() - 1
    threads_updated = max(1, min(threads, available_threads))
    return threads_updated


def generate_report(name: str) -> None:
    report_generator.report(name)
    print(
        f"\N{party popper} Finished! Open {config.DAJIN_RESULTS_DIR}/{name} to see the report.",
        file=sys.stderr,
    )


################################################################################
# Single mode
################################################################################


def execute_single_mode(arguments: dict[str]):
    input_validator.validate_files(arguments["sample"], arguments["control"], arguments["allele"])
    if "genome" in arguments:
        arguments.update(input_validator.validate_genome_and_fetch_urls(arguments["genome"]))
    core.execute_control(arguments)
    core.execute_sample(arguments)
    generate_report(arguments["name"])


################################################################################
# Batch mode
################################################################################


def validate_columns(columns: list, filepath: str) -> None:
    """Validate the columns of a batch file."""
    required_columns = ["sample", "control", "allele", "name"]
    accepted_columns = ["sample", "control", "allele", "name", "genome"]

    if not set(required_columns).issubset(set(columns)):
        raise ValueError(f"{filepath} must contain {', '.join(required_columns)} in the header")

    if not set(columns).issubset(accepted_columns):
        raise ValueError(f"Accepted header names of {filepath} are {', '.join(accepted_columns)}.")


def _batched(iterable, chunk_size: int) -> Generator(tuple):
    iterator = iter(iterable)
    while True:
        chunk = tuple(islice(iterator, chunk_size))
        if not chunk:
            return
        yield chunk


def run_multiprocess(function, arguments: list[dict], num_workers: int = 1) -> None:
    logger = logging.getLogger(__name__)

    q = Queue()

    def target(arg):
        try:
            function(arg)
        except Exception:
            q.put(traceback.format_exc())
            sys.exit(1)

    arguments_batched = _batched(arguments, num_workers)
    for args in arguments_batched:
        processes = [Process(target=target, args=(arg,)) for arg in args]

        for p in processes:
            p.start()

        for p in processes:
            p.join()
            if p.exitcode == 1:
                logger.error(f"An unexpected error occurred: {q.get()}")
                sys.exit(1)


def execute_batch_mode(arguments: dict[str]):
    path_batchfile = io.convert_to_posix(arguments["file"])

    if not Path(path_batchfile).exists():
        raise FileNotFoundError(f"'{path_batchfile}' does not exist.")

    inputs = io.load_file(path_batchfile)

    # Validate Column of the batch file
    columns = inputs[0]
    validate_columns(columns, path_batchfile)

    # Validate contents and fetch genome urls
    contents = inputs[1:]
    cache_urls_genome = dict()
    index_of_name = columns.index("name")
    contents.sort(key=lambda x: x[index_of_name])
    for _, groups in groupby(contents, key=lambda x: x[index_of_name]):
        for group in groups:
            args = {h: g for h, g in zip(columns, group)}
            # validate contents in the batch file
            input_validator.validate_files(args["sample"], args["control"], args["allele"])
            # validate genome and fetch urls
            if "genome" in args and args["genome"] not in cache_urls_genome:
                urls_genome = input_validator.validate_genome_and_fetch_urls(args["genome"])
                cache_urls_genome[args["genome"]] = urls_genome

    ##############################################################################
    # Perform DAJIN2
    ##############################################################################

    index_of_name = columns.index("name")
    contents.sort(key=lambda x: x[index_of_name])
    contents_sample = []
    for name, groups in groupby(contents, key=lambda x: x[index_of_name]):
        groups = list(groups)

        # Handle controls
        contents_control = []
        for group in groups:
            args = {h: g for h, g in zip(columns, group)}
            args["threads"] = arguments["threads"]
            args["sample"] = args["control"]
            if "genome" in args:
                args.update(cache_urls_genome[args["genome"]])
            contents_control.append(args)
        contents_control_unique = [dict(item) for item in set(frozenset(d.items()) for d in contents_control)]
        run_multiprocess(core.execute_control, contents_control_unique, arguments["threads"])

        # Handle samples
        contents_sample = []
        for group in groups:
            args = {h: g for h, g in zip(columns, group)}
            args["threads"] = arguments["threads"]
            if args["sample"] == args["control"]:
                continue
            if "genome" in args:
                args.update(cache_urls_genome[args["genome"]])
            contents_sample.append(args)
        contents_sample_unique = [dict(item) for item in set(frozenset(d.items()) for d in contents_sample)]
        run_multiprocess(core.execute_sample, contents_sample_unique, arguments["threads"])
        # Finish
        generate_report(name)


def execute():
    parser = argparse.ArgumentParser()

    ###############################################################################
    # Single mode
    ###############################################################################

    parser.add_argument("-s", "--sample", type=str, help="Full path to a sample FASTQ file")
    parser.add_argument("-c", "--control", type=str, help="Full path to a control FASTQ file")
    parser.add_argument("-a", "--allele", type=str, help="Full path to a FASTA file")
    parser.add_argument("-n", "--name", type=str, help="Output directory name", default="DAJIN2-results")
    parser.add_argument(
        "-g", "--genome", type=str, default="", help="Reference genome ID (e.g hg38, mm39) [default: '']"
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
        arguments["threads"] = update_threads(int(args.threads))
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
            raise ValueError("the following arguments are required: -s/--sample")
        if args.control is None:
            raise ValueError("the following arguments are required: -c/--control")
        if args.allele is None:
            raise ValueError("the following arguments are required: -a/--allele")
        if args.name is None:
            raise ValueError("the following arguments are required: -n/--name")
        arguments = dict()
        arguments["sample"] = args.sample
        arguments["control"] = args.control
        arguments["allele"] = args.allele
        arguments["name"] = args.name
        arguments["genome"] = args.genome
        arguments["threads"] = update_threads(int(args.threads))
        arguments["debug"] = args.debug

        logger = logging.getLogger(__name__)
        try:
            execute_single_mode(arguments)
        except Exception as e:
            logger.error(f"An unexpected error occurred: {e}")
            sys.exit(1)


if __name__ == "__main__":
    execute()
