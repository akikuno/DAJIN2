import argparser
import shutil
import sys
import os
import re
from typing import List


def check_deps(dependencies: List[str]) -> None:
    for dep in dependencies:
        if not shutil.which(dep):
            print(f"{dep} is not found. Please install it.")
            sys.exit(1)


def make_dirs(main: str, subfolders: List[str]) -> None:
    for sub in subfolders:
        dir = os.path.join(main, sub)
        os.makedirs(dir, exist_ok=True)


def main():
    check_deps(["minimap2", "samtools"])

    make_dirs(".tmpDAJIN", ["fasta", "sam"])

    sample, control, allele, output, genome, debug, threads = argparser.parse()

    split_fasta()
    if debug is False:
        shutil.rmtree(".tmpDAJIN")


if __name__ == "__main__":
    main()
