#!/usr/bin/env python
import argparse
import os

description = "DAJIN: Genotyping software using Nanopore long-read sequencer for genome-edited samples."`

default_proc = max(len(os.sched_getaffinity(0))-1, 1)

parser = argparse.ArgumentParser(description=description)

parser.add_argument("-a", "--alleles", required=True)
parser.add_argument("-c", "--control", required=True)
parser.add_argument("-s", "--sample", required=True)
parser.add_argument("-s", "--sample", required=True)
parser.add_argument("-g", "--genome")
parser.add_argument("-o", "--output", default="DAJIN_results")
parser.add_argument("-t", "--threads", default=default_proc)
