from __future__ import annotations

import os

os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"

import argparse
import importlib.metadata
import logging
import shutil
import sys
from copy import deepcopy
from itertools import groupby
from pathlib import Path

from DAJIN2 import gui, view
from DAJIN2.core import core
from DAJIN2.utils import config, input_validator, io, multiprocess, report_generator

try:
    DAJIN_VERSION = importlib.metadata.version("DAJIN2")
except importlib.metadata.PackageNotFoundError:
    DAJIN_VERSION = "counld not be retrieved"


def generate_report(name: str, logger: logging.Logger) -> None:
    report_generator.report(name)
    logger.info(f"\N{PARTY POPPER} Finished! Open {config.DAJIN_RESULTS_DIR}/{name} to see the report.")


################################################################################
# Single mode
################################################################################


def execute_single_mode(arguments: dict[str]):
    # Set logging to export log to stderr and file
    path_logfile = config.get_logfile()
    logger = config.set_logging(path_logfile)
    logger.info(f"\N{RUNNER} Start running DAJIN2 version {DAJIN_VERSION}")
    logger.info(f"\N{PERSONAL COMPUTER} {' '.join(sys.argv)}")

    # Validate input files
    input_validator.validate_files(arguments["sample"], arguments["control"], arguments["allele"])
    if arguments.get("genome"):
        arguments.update(input_validator.validate_genome_and_fetch_urls(arguments["genome"]))

    # Run DAJIN2
    core.execute_control(arguments)
    core.execute_sample(arguments)

    # Finish call
    generate_report(arguments["name"], logger)
    shutil.move(path_logfile, Path("DAJIN_Results", arguments["name"]))

    # Remove temporary files
    if not arguments["debug"]:
        shutil.rmtree(Path("DAJIN_Results", ".tempdir", arguments["name"]))


################################################################################
# Batch mode
################################################################################


def validate_headers_of_batch_file(headers: set[str], filepath: str) -> None:
    """Validate the headers of a batch file."""
    required_headers = {"sample", "control", "allele", "name"}
    accepted_headers = {"sample", "control", "allele", "name", "genome"}

    if not required_headers.issubset(headers):
        raise ValueError(f'{filepath} must contain "sample", "control", "allele" and "name" in the header')

    if not headers.issubset(accepted_headers):
        raise ValueError(
            f'Accepted header names of {filepath} are "sample", "control", "allele", "name", or "genome".'
        )


def create_argument_dict(args: dict, cache_urls_genome: dict, is_control: bool) -> dict[str, str]:
    """Create a dictionary of arguments from the given headers and group."""
    args_update = deepcopy(args)

    args_update["threads"] = 1  # Set the number of threads to 1 for batch mode

    # Assign the "sample" field depending on whether it's a control or not
    if is_control:
        args_update["sample"] = args_update["control"]
    else:
        if args_update["sample"] == args_update["control"]:
            return {}  # Return an empty dict to indicate a skipped group

    if args_update.get("genome"):
        args_update.update(cache_urls_genome[args_update["genome"]])

    return args_update


def run_DAJIN2(
    groups: list[dict[str, str]], cache_urls_genome: dict, is_control: bool = True, num_workers: int = 1
) -> None:
    contents = []
    for args in groups:
        args = create_argument_dict(args, cache_urls_genome, is_control)
        if args:  # Add args to contents only if it's not an empty dict
            contents.append(args)

    # Return a list of unique dictionaries
    contents_unique = [dict(item) for item in {frozenset(d.items()) for d in contents}]

    contents_unique.sort(key=lambda x: x["sample"])

    if is_control:
        multiprocess.run(core.execute_control, contents_unique, num_workers)
    else:
        multiprocess.run(core.execute_sample, contents_unique, num_workers)


def execute_batch_mode(arguments: dict[str]):
    path_batchfile = io.convert_to_posix(arguments["file"])

    if not Path(path_batchfile).exists():
        raise FileNotFoundError(f"'{path_batchfile}' does not exist.")

    records = io.load_batchfile(path_batchfile)

    # Validate Column of the batch file
    headers = set(records[0].keys())
    validate_headers_of_batch_file(headers, path_batchfile)

    # Validate contents and fetch genome urls
    cache_urls_genome = {}
    records.sort(key=lambda x: x["name"])
    for _, groups in groupby(records, key=lambda x: x["name"]):
        for args in groups:
            # Validate contents in the batch file
            input_validator.validate_files(args["sample"], args["control"], args["allele"])
            # Validate genome and fetch urls
            if args.get("genome") and args["genome"] not in cache_urls_genome:
                urls_genome = input_validator.validate_genome_and_fetch_urls(args["genome"])
                cache_urls_genome[args["genome"]] = urls_genome

    # Run DAJIN2
    for name, groups in groupby(records, key=lambda x: x["name"]):
        groups: list[dict[str, str]] = list(groups)
        # Set logging to export log to stderr and file
        config.reset_logging()
        path_logfile = config.get_logfile()
        logger = config.set_logging(path_logfile)

        logger.info(f"\N{RUNNER} Start running DAJIN2 version {DAJIN_VERSION}")
        logger.info(f"\N{PERSONAL COMPUTER} {' '.join(sys.argv)}")
        logger.info(f"\N{LEFT-POINTING MAGNIFYING GLASS} Handling {name}")

        # Run DAJIN2
        run_DAJIN2(groups, cache_urls_genome, is_control=True, num_workers=arguments["threads"])
        run_DAJIN2(groups, cache_urls_genome, is_control=False, num_workers=arguments["threads"])

        # Finish call
        generate_report(name, logger)
        shutil.move(path_logfile, Path("DAJIN_Results", name))
        if not arguments["debug"]:
            shutil.rmtree(Path("DAJIN_Results", ".tempdir", name))


def execute():
    parser = argparse.ArgumentParser()

    ###############################################################################
    # Single mode
    ###############################################################################

    parser.add_argument("-s", "--sample", type=str, help="Path to a sample directory including FASTQ file")
    parser.add_argument("-c", "--control", type=str, help="Path to a control directory including FASTQ file")
    parser.add_argument("-a", "--allele", type=str, help="Path to a FASTA file")
    parser.add_argument("-n", "--name", type=str, help="Output directory name", default="Results")
    parser.add_argument(
        "-g", "--genome", type=str, default="", help="Reference genome ID (e.g hg38, mm39) [default: '']"
    )
    parser.add_argument("-t", "--threads", type=int, default=1, help="Number of threads [default: 1]")
    parser.add_argument("-v", "--version", action="version", version=f"DAJIN2 version {DAJIN_VERSION}")
    parser.add_argument("-d", "--debug", action="store_true", help=argparse.SUPPRESS)

    ###############################################################################
    # Batch mode
    ###############################################################################

    def batchmode(args):
        arguments = {}
        arguments["file"] = args.file
        arguments["threads"] = input_validator.update_threads(int(args.threads))
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
        arguments = {}
        arguments["sample"] = args.sample
        arguments["control"] = args.control
        arguments["allele"] = args.allele
        arguments["name"] = args.name
        arguments["genome"] = args.genome
        arguments["threads"] = input_validator.update_threads(int(args.threads))
        arguments["debug"] = args.debug

        execute_single_mode(arguments)


if __name__ == "__main__":
    execute()
