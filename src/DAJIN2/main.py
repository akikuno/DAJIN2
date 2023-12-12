from __future__ import annotations

from DAJIN2.utils import config

# prevent BLAS from using all cores
config.set_single_threaded_blas()

import sys
import shutil
import argparse
from pathlib import Path
from itertools import groupby

from DAJIN2 import gui, view
from DAJIN2.core import core
from DAJIN2.utils import config, io, report_generator, input_validator, multiprocess


DAJIN_VERSION = "0.3.4"


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
    # set logging to export log to stderr and file
    path_logfile = config.get_logfile()
    config.set_logging(path_logfile)
    input_validator.validate_files(arguments["sample"], arguments["control"], arguments["allele"])
    if arguments.get("genome"):
        arguments.update(input_validator.validate_genome_and_fetch_urls(arguments["genome"]))
    # run DAJIN2
    core.execute_control(arguments)
    core.execute_sample(arguments)
    # finish
    generate_report(arguments["name"])
    shutil.move(path_logfile, Path("DAJIN_Results", arguments["name"]))
    if not arguments["debug"]:
        shutil.rmtree(Path("DAJIN_Results", ".tempdir", arguments["name"]))


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


def create_argument_dict(columns: list, group: list, cache_urls_genome: dict, is_control: bool) -> dict:
    """Create a dictionary of arguments from the given columns and group."""
    args = dict(zip(columns, group))
    args["threads"] = 1  # Set the number of threads to 1 for batch mode

    # Assign the "sample" field depending on whether it's a control or not
    if is_control:
        args["sample"] = args["control"]
    else:
        if args["sample"] == args["control"]:
            return {}  # Return an empty dict to indicate a skipped group

    if args.get("genome"):
        args.update(cache_urls_genome[args["genome"]])

    return args


def run_DAJIN2(
    groups: list, columns: list, cache_urls_genome: dict, is_control: bool = True, num_workers: int = 1
) -> None:
    contents = []
    for group in groups:
        args = create_argument_dict(columns, group, cache_urls_genome, is_control)
        if args:  # Add args to contents only if it's not an empty dict
            contents.append(args)

    # Return a list of unique dictionaries
    contents_unique = [dict(item) for item in set(frozenset(d.items()) for d in contents)]

    contents_unique.sort(key=lambda x: x["sample"])

    if is_control:
        multiprocess.run(core.execute_control, contents_unique, num_workers)
    else:
        multiprocess.run(core.execute_sample, contents_unique, num_workers)


def execute_batch_mode(arguments: dict[str]):
    path_batchfile = io.convert_to_posix(arguments["file"])

    if not Path(path_batchfile).exists():
        raise FileNotFoundError(f"'{path_batchfile}' does not exist.")

    inputs = io.load_batchfile(path_batchfile)

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
            args = dict(zip(columns, group))
            # validate contents in the batch file
            input_validator.validate_files(args["sample"], args["control"], args["allele"])
            # validate genome and fetch urls
            if args.get("genome") and args["genome"] not in cache_urls_genome:
                urls_genome = input_validator.validate_genome_and_fetch_urls(args["genome"])
                cache_urls_genome[args["genome"]] = urls_genome
    for name, groups in groupby(contents, key=lambda x: x[index_of_name]):
        # set logging to export log to stderr and file
        path_logfile = config.get_logfile()
        config.set_logging(path_logfile)
        groups = list(groups)
        # Run DAJIN2
        run_DAJIN2(groups, columns, cache_urls_genome, is_control=True, num_workers=arguments["threads"])
        run_DAJIN2(groups, columns, cache_urls_genome, is_control=False, num_workers=arguments["threads"])
        # Finish
        generate_report(name)
        shutil.move(path_logfile, Path("DAJIN_Results", name))
        if not arguments["debug"]:
            shutil.rmtree(Path("DAJIN_Results", ".tempdir", name))


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
    parser.add_argument("-v", "--version", action="version", version=f"DAJIN2 version {DAJIN_VERSION}")
    parser.add_argument("-d", "--debug", action="store_true", help=argparse.SUPPRESS)

    ###############################################################################
    # Batch mode
    ###############################################################################

    def batchmode(args):
        arguments = dict()
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
        arguments = dict()
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
