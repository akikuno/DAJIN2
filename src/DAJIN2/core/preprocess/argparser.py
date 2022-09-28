import argparse
import os


def parse():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-s", "--sample", required=True, metavar="<sample.fastq>", help="Give the full path to a control FASTQ file",
    )
    parser.add_argument(
        "-c", "--control", required=True, metavar="<control.fastq>", help="Give the full path to a control FASTQ file",
    )
    parser.add_argument(
        "-a", "--allele", required=True, metavar="<allele.fasta>", help="Give the full path to allele FASTA file",
    )
    parser.add_argument(
        "-o",
        "--output",
        metavar="<output directory>",
        default="DAJIN_results",
        type=str,
        help="Name of the output directory",
    )
    parser.add_argument(
        "-g", "--genome", default="", metavar="<genome id>", help="Reference genome ID (e.g hg38, mm10)",
    )
    parser.add_argument(
        "-t", "--threads", default=1, type=int, metavar="<threads>", help="Number of threads [default: 1]",
    )
    parser.add_argument("-d", "--debug", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("-v", "--version", action="version", version="DAJIN2 version 1.0.0")
    args = parser.parse_args()

    os_cpus = int(os.cpu_count())
    if args.threads > os_cpus:
        threads = os_cpus
    elif args.threads < 1:
        threads = 1
    else:
        threads = args.threads

    return (
        args.sample,
        args.control,
        args.allele,
        args.output,
        args.genome,
        args.debug,
        threads,
    )
