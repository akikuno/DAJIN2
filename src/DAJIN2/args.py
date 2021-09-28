import argparse
import os


def parse():
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--sample",
                        required=True,
                        metavar="<sample.fastq>",
                        help="Give the full path to a control FASTQ file")
    parser.add_argument("-c", "--control",
                        required=True,
                        metavar="<control.fasta>",
                        help="Give the full path to a control FASTA file")
    parser.add_argument("-o", "--output",
                        metavar="<output directory>",
                        default="DAJIN_results",
                        type=str,
                        help="Name of the output directory")
    parser.add_argument("-g", "--genome",
                        metavar="<genome id>",
                        help="Reference genome ID (e.g hg38, mm10)")
    parser.add_argument("-d", "--debug",
                        action="store_true",
                        help="Retain all intermediate files")
    parser.add_argument("-t", "--threads",
                        default=1,
                        type=int,
                        metavar="<threads>",
                        help="Number of threads [default: 1]",)
    parser.add_argument('-v', '--version',
                        action='version',
                        version='DAJIN version 2.0.0')
    args = parser.parse_args()

    os.makedirs(args.output, exist_ok=True)
    os_cpus = int(os.cpu_count())  # len(os.sched_getaffinity(0))
    if args.threads > os_cpus:
        threads = os_cpus
    elif args.threads < 1:
        threads = 1
    else:
        threads = args.threads
    return (args.sample, args.control, args.output,
            args.genome, args.debug, threads)
