from __future__ import annotations

import os
import sys
import logging
import argparse

from pathlib import Path
from itertools import groupby
from DAJIN2 import gui, view
from DAJIN2.core import core
from DAJIN2.utils import config, io, report_generator, input_validator, multiprocess


VERSION = "0.3.2"

# prevent BLAS from using all cores
config.set_single_threaded_blas()

# set logging to export log to stderr and file
config.set_logging()


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


def validate_columns_of_batch_file(columns: list, filepath: str) -> None:
    """Validate the columns of a batch file."""
    required_columns = ["sample", "control", "allele", "name"]
    accepted_columns = ["sample", "control", "allele", "name", "genome"]

    if not set(required_columns).issubset(set(columns)):
        raise ValueError(f"{filepath} must contain {', '.join(required_columns)} in the header")

    if not set(columns).issubset(accepted_columns):
        raise ValueError(f"Accepted header names of {filepath} are {', '.join(accepted_columns)}.")


def execute_batch_mode(arguments: dict[str]):
    path_batchfile = io.convert_to_posix(arguments["file"])

    if not Path(path_batchfile).exists():
        raise FileNotFoundError(f"'{path_batchfile}' does not exist.")

    inputs = io.load_file(path_batchfile)

    # Validate Column of the batch file
    columns = inputs[0]
    validate_columns_of_batch_file(columns, path_batchfile)

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

    def run_DAJIN2(groups: list, is_control: bool = True, num_workers: int = 1) -> None:
        contents = []
        for group in groups:
            args = {h: g for h, g in zip(columns, group)}
            # Set the number of threads to 1 for batch mode
            args["threads"] = 1
            # Assign the "sample" field depending on whether it's a control or not
            if is_control:
                args["sample"] = args["control"]
            else:
                if args["sample"] == args["control"]:
                    continue
            if "genome" in args:
                args.update(cache_urls_genome[args["genome"]])
            contents.append(args)
        # return a list of unique dictionaries from a list of dictionaries
        contents_unique = [dict(item) for item in set(frozenset(d.items()) for d in contents)]
        if is_control:
            multiprocess.run(core.execute_control, contents_unique, num_workers)
        else:
            multiprocess.run(core.execute_sample, contents_unique, num_workers)

    for name, groups in groupby(contents, key=lambda x: x[index_of_name]):
        groups = list(groups)
        # Run DAJIN2
        run_DAJIN2(groups, is_control=True, num_workers=arguments["threads"])
        run_DAJIN2(groups, is_control=False, num_workers=arguments["threads"])
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
