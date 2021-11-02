import argparser
import shutil
import sys
import os
import tempfile

# Custom modules
import cache
from preprocess import mids
from preprocess import mapping


def check_dependencies(dependencies: list) -> None:
    for dep in dependencies:
        if not shutil.which(dep):
            print(f"{dep} is not found. Please install it.")
            sys.exit(1)


def make_directories(maindir: str, subdirs: list) -> None:
    for sub in subdirs:
        dir = os.path.join(maindir, sub)
        os.makedirs(dir, exist_ok=True)


def main():
    ########################################################################
    # Setting
    ########################################################################
    check_dependencies(["minimap2", "samtools"])

    TMPDIR = ".tmpDAJIN"
    os.makedirs(TMPDIR, exist_ok=True)
    make_directories(TMPDIR, ["fasta", "sam"])
    TMPDIR_PATHS = {_: os.path.join(TMPDIR, _) for _ in os.listdir(TMPDIR)}

    CACHEDIR = os.path.join(tempfile.gettempdir(), "DAJIN")
    os.makedirs(CACHEDIR, exist_ok=True)

    ########################################################################
    # Argparse
    ########################################################################

    sample, control, allele, output, genome, debug, threads = argparser.parse()


    if cache.exists(control, CACHEDIR) and cache.is_same_filesize(control, CACHEDIR):
        IS_CACHED = True
    else:
        shutil.copy(control, CACHEDIR)
        IS_CACHED = False

    ########################################################################
    # minimap2
    ########################################################################

    mapping.split_fasta(allele, TMPDIR_PATHS["fasta"])

    if IS_CACHED:
        cache.load(CACHEDIR, TMPDIR_PATHS["sam"])
    else:
        mapping.minimap2(control, TMPDIR_PATHS["fasta"], TMPDIR_PATHS["sam"], threads)
        cache.save(control, TMPDIR_PATHS["sam"], CACHEDIR)

    mapping.minimap2(sample, TMPDIR_PATHS["fasta"], TMPDIR_PATHS["sam"], threads)

    ########################################################################
    # MIDS conversion
    ########################################################################
    
    SAMDIR = os.path.join("tests", "samTomids", "input")
    SAMFILES = [os.path.join(SAMDIR, _) for _ in os.listdir(SAMDIR)]

    samfile = SAMFILES[-3]

    mids_csv = mids.samfile_to_mids(samfile)
    print(mids_csv)

    if debug is False:
        shutil.rmtree(".tmpDAJIN")


if __name__ == "__main__":
    main()
