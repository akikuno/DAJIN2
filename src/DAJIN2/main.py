import argparser
import shutil
import sys
import os
import mapping
from typing import List
import tempfile
import glob


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

    TMPDIR = ".tmpDAJIN"
    make_dirs(TMPDIR, ["fasta", "sam"])
    TMPDIRS = {_: os.path.join(TMPDIR, _) for _ in os.listdir(TMPDIR)}

    CACHEDIR = os.path.join(tempfile.gettempdir(), "DAJIN")
    os.makedirs(CACHEDIR, exist_ok=True)

    sample, control, allele, output, genome, debug, threads = argparser.parse()

    mapping.split_fasta(allele, TMPDIRS["fasta"])
    mapping.minimap2(control, TMPDIRS["fasta"], TMPDIRS["sam"], threads)
    mapping.minimap2(sample, TMPDIRS["fasta"], TMPDIRS["sam"], threads)
    if debug is False:
        shutil.rmtree(".tmpDAJIN")


if __name__ == "__main__":
    main()
