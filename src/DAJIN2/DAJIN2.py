import os
import argparse

from .singlemode import single
from .batchmode import batch


def main():
    parser = argparse.ArgumentParser()

    ###############################################################################
    # Single mode
    ###############################################################################

    parser.add_argument("-s", "--sample", type=str, help="Give the full path to a control FASTQ file")
    parser.add_argument("-c", "--control", type=str, help="Give the full path to a control FASTQ file")
    parser.add_argument("-a", "--allele", type=str, help="Give the full path to allele FASTA file")
    parser.add_argument("-n", "--name", type=str, help="Output name of the report")
    parser.add_argument(
        "-g", "--genome", type=str, default="", help="Reference genome ID (e.g hg38, mm10) [default: '']"
    )
    parser.add_argument("-t", "--threads", type=int, default=1, help="Number of threads [default: 1]")
    parser.add_argument("-v", "--version", action="version", version="DAJIN2 version 1.0.0")
    parser.add_argument("-d", "--debug", action="store_true", help=argparse.SUPPRESS)

    ###############################################################################
    # Batch mode
    ###############################################################################

    subparser = parser.add_subparsers()
    parser_batch = subparser.add_parser("batch", help="DAIJN2 batch mode")
    parser_batch.add_argument("-f", "--file", required=True, type=str, help="batch.csv")
    parser_batch.add_argument("-t", "--threads", default=1, type=int, help="Number of threads [default: 1]")
    parser_batch.add_argument("-d", "--debug", action="store_true", help=argparse.SUPPRESS)

    args = parser.parse_args()

    ###############################################################################
    # Parse arguments
    ###############################################################################

    # Update threads
    os_cpus = int(os.cpu_count())
    threads = max(1, args.threads)
    if args.threads > os_cpus:
        threads = max(1, os_cpus - 1)

    try:
        _ = args.file
    except AttributeError:
        if args.sample is None:
            raise AttributeError("-s/--sample is required")
        if args.control is None:
            raise AttributeError("-c/--control is required")
        if args.allele is None:
            raise AttributeError("-a/--allele is required")
        if args.name is None:
            raise AttributeError("-n/--name is required")
        arguments = dict()
        arguments["sample"] = args.sample
        arguments["control"] = args.control
        arguments["allele"] = args.allele
        arguments["name"] = args.name
        arguments["genome"] = args.genome
        arguments["threads"] = threads
        arguments["debug"] = args.debug
        single.single(arguments)
    else:
        arguments = dict()
        arguments["file"] = args.file
        arguments["threads"] = threads
        arguments["debug"] = args.debug
        batch.batch(arguments)


if __name__ == "__main__":
    main()
