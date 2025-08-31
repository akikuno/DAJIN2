from __future__ import annotations

import os

os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"

import argparse
import importlib.metadata
import json
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


def cache_bed_coordinates(name: str, genome_coordinates: dict, logger: logging.Logger) -> None:
    """Cache BED file coordinates to genome_coordinates.jsonl."""
    try:
        # Create cache directory structure: .tempdir/{NAME}/cache/
        cache_dir = Path(config.DAJIN_RESULTS_DIR) / ".tempdir" / name / "cache"
        cache_dir.mkdir(parents=True, exist_ok=True)

        # Save to genome_coordinates.jsonl
        cache_file = cache_dir / "genome_coordinates.jsonl"
        with open(cache_file, "w") as f:
            json.dump(genome_coordinates, f)
            f.write("\n")

        logger.debug(f"Cached BED coordinates to {cache_file}")

    except Exception as e:
        logger.warning(f"Failed to cache BED coordinates: {e}")


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

    # Validate genome coordinates (BED file takes precedence over genome ID)
    if arguments["genome_coordinate"]:
        input_validator.validate_bed_file(arguments["genome_coordinate"])
    elif arguments["genome"]:
        # genome ID provided - use GGGenome and UCSC lookup
        # Here, two processes are executed exceptionally: validating the server and obtaining valid server information.
        # The reason is that splitting these processes would require accessing the server twice, which would increase the server load.
        arguments |= input_validator.get_available_servers()

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
    accepted_headers = {"sample", "control", "allele", "name", "genome", "genome_coordinate", "bed"}

    if not required_headers.issubset(headers):
        raise ValueError(
            f'{filepath} must contain "sample", "control", "allele" and "name" in the header. Currently, it contains {headers}.'
        )

    if not headers.issubset(accepted_headers):
        raise ValueError(
            f'Accepted header names of {filepath} are "sample", "control", "allele", "name", "genome", "genome_coordinate", or "bed". Currently, it contains {headers}.'
        )


def create_argument_dict(args: dict, genome_with_server_urls: dict[str, str], is_control: bool) -> dict[str, str]:
    """Create a dictionary of arguments from the given headers and group."""
    args_update = deepcopy(args)

    args_update["threads"] = 1  # Set the number of threads to 1 for batch mode

    # Assign the "sample" field depending on whether it's a control or not
    if is_control:
        args_update["sample"] = args_update["control"]
    else:
        if args_update["sample"] == args_update["control"]:
            return {}  # Return an empty dict to indicate a skipped group

    # Handle genome coordinates (BED file takes precedence over genome ID)
    bed_file = args_update.get("genome_coordinate") or args_update.get("bed")

    if bed_file:
        # Ensure genome_coordinate is set for consistency
        args_update["genome_coordinate"] = bed_file
    elif args_update.get("genome"):
        args_update |= genome_with_server_urls

    return args_update


def run_DAJIN2(
    groups: list[dict[str, str]],
    genome_with_server_urls: dict[str, str],
    is_control: bool = True,
    num_workers: int = 1,
    no_filter: bool = False,
) -> None:
    contents = []
    for args in groups:
        args = create_argument_dict(args, genome_with_server_urls, is_control)
        if args:  # Add args to contents only if it's not an empty dict
            args["no_filter"] = no_filter  # Add no_filter to each args dict
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
    genome_with_server_urls = {}
    records.sort(key=lambda x: x["name"])
    for _, groups in groupby(records, key=lambda x: x["name"]):
        for args in groups:
            # Validate contents in the batch file
            input_validator.validate_files(args["sample"], args["control"], args["allele"])

            # Validate BED file if provided
            if args.get("genome_coordinate"):
                input_validator.validate_file_existence(args["genome_coordinate"])

            # Validate genome and fetch urls
            if args.get("genome") and args["genome"] not in genome_with_server_urls:
                genome_with_server_urls["genome"] = args["genome"]
                genome_with_server_urls |= input_validator.get_available_servers()

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
        run_DAJIN2(
            groups,
            genome_with_server_urls,
            is_control=True,
            num_workers=arguments["threads"],
            no_filter=arguments["no_filter"],
        )
        run_DAJIN2(
            groups,
            genome_with_server_urls,
            is_control=False,
            num_workers=arguments["threads"],
            no_filter=arguments["no_filter"],
        )

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
        "-g", "--genome", type=str, default=None, help="Reference genome ID (e.g hg38, mm39) [default: '']"
    )
    parser.add_argument(
        "-b",
        "--bed",
        type=str,
        default=None,
        dest="genome_coordinate",
        help="Path to BED6 file containing genomic coordinates [default: '']",
    )
    parser.add_argument("-t", "--threads", type=int, default=1, help="Number of threads [default: 1]")
    parser.add_argument(
        "--no-filter", action="store_true", help="Disable minor allele filtering (keep alleles <0.5%%)"
    )
    parser.add_argument("-v", "--version", action="version", version=f"DAJIN2 version {DAJIN_VERSION}")
    parser.add_argument("-d", "--debug", action="store_true", help=argparse.SUPPRESS)

    ###############################################################################
    # Batch mode
    ###############################################################################

    def batchmode(args):
        arguments = {}
        arguments["file"] = args.file
        arguments["threads"] = input_validator.update_threads(int(args.threads))
        arguments["no_filter"] = args.no_filter
        arguments["debug"] = args.debug
        execute_batch_mode(arguments)

    subparser = parser.add_subparsers()
    parser_batch = subparser.add_parser("batch", help="DAIJN2 batch mode")
    parser_batch.add_argument("-f", "--file", required=True, type=str, help="CSV or Excel file.")
    parser_batch.add_argument("-t", "--threads", default=1, type=int, help="Number of threads [default: 1]")
    parser_batch.add_argument(
        "--no-filter", action="store_true", help="Disable minor allele filtering (keep alleles <0.5%%)"
    )
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
            parser.error("the following arguments are required: -s/--sample")
        if args.control is None:
            parser.error("the following arguments are required: -c/--control")
        if args.allele is None:
            parser.error("the following arguments are required: -a/--allele")
        if args.name is None:
            parser.error("the following arguments are required: -n/--name")

        arguments = {}
        arguments["sample"] = args.sample
        arguments["control"] = args.control
        arguments["allele"] = args.allele
        arguments["name"] = args.name
        arguments["genome"] = args.genome
        arguments["genome_coordinate"] = args.genome_coordinate
        arguments["threads"] = input_validator.update_threads(int(args.threads))
        arguments["no_filter"] = args.no_filter
        arguments["debug"] = args.debug

        execute_single_mode(arguments)


if __name__ == "__main__":
    execute()
