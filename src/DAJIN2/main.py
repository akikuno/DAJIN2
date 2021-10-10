import argparser
import shutil
import sys
import os
import tempfile

# Custom modules
import cache
from preprocess import mids
from preprocess import mapping


def check_deps(dependencies: list) -> None:
    for dep in dependencies:
        if not shutil.which(dep):
            print(f"{dep} is not found. Please install it.")
            sys.exit(1)


def make_dirs(maindir: str, subdirs: list) -> None:
    for sub in subdirs:
        dir = os.path.join(maindir, sub)
        os.makedirs(dir, exist_ok=True)


def main():
    check_deps(["minimap2", "samtools"])

    TMPDIR = ".tmpDAJIN"
    os.makedirs(TMPDIR, exist_ok=True)
    make_dirs(TMPDIR, ["fasta", "sam"])
    TMPDIRS = {_: os.path.join(TMPDIR, _) for _ in os.listdir(TMPDIR)}

    CACHEDIR = os.path.join(tempfile.gettempdir(), "DAJIN")
    os.makedirs(CACHEDIR, exist_ok=True)

    sample, control, allele, output, genome, debug, threads = argparser.parse()

    if cache.exists(control, CACHEDIR) and cache.is_same_size(control, CACHEDIR):
        IS_CACHED = True
    else:
        shutil.copy(control, CACHEDIR)
        IS_CACHED = False

    ########################################################################
    # minimap2
    ########################################################################

    mapping.split_fasta(allele, TMPDIRS["fasta"])

    if IS_CACHED:
        cache.load(CACHEDIR, TMPDIRS["sam"])
    else:
        mapping.minimap2(control, TMPDIRS["fasta"], TMPDIRS["sam"], threads)
        cache.save(control, TMPDIRS["sam"], CACHEDIR)

    mapping.minimap2(sample, TMPDIRS["fasta"], TMPDIRS["sam"], threads)

    ########################################################################
    # MIDS conversion
    ########################################################################
    SAMDIR = os.path.join("tests", "samTomids", "input")
    SAMFILES = [os.path.join(SAMDIR, _) for _ in os.listdir(SAMDIR)]

    samfile = SAMFILES[-3]

    mids = mids.samfile_to_mids(samfile)
    print(mids)

    if debug is False:
        shutil.rmtree(".tmpDAJIN")


if __name__ == "__main__":
    main()
