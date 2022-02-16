import shutil
import sys
import os
import tempfile

# Custom modules
from src.DAJIN2.utils import cache_control
from src.DAJIN2.utils import argparser
from src.DAJIN2.preprocess import mapping
from src.DAJIN2.preprocess import midsconv
from src.DAJIN2.preprocess import midsqc

# For development
# import importlib


#! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# def main():
########################################################################
# Setting
########################################################################


def check_dependencies(dependencies: list) -> None:
    for dep in dependencies:
        if not shutil.which(dep):
            print(f"{dep} is not found. Please install it.")
            sys.exit(1)


check_dependencies(["minimap2", "samtools"])


def make_directories(maindir: str, subdirs: list) -> None:
    os.makedirs(maindir, exist_ok=True)
    for sub in subdirs:
        dir = os.path.join(maindir, sub)
        os.makedirs(dir, exist_ok=True)


TMPDIR = ".tmpDAJIN"
SUBDIRS = ["fasta", "sam", "midsconv", "midsqc"]
for sub in SUBDIRS:
    dir = os.path.join(TMPDIR, sub)
    os.makedirs(dir, exist_ok=True)

TMPDIR_PATHS = {
    dirname: os.path.join(TMPDIR, dirname) for dirname in os.listdir(TMPDIR)
}

########################################################################
# Argument parse
########################################################################

sample, control, allele, output, genome, debug, threads = argparser.parse()

# sample, control, allele, output, genome, debug, threads = (
#     "examples/pm-tyr/barcode31.fq.gz",
#     "examples/pm-tyr/barcode32.fq.gz",
#     "examples/pm-tyr/design_tyr.fa",
#     "DAJIN_results",
#     "mm10",
#     True,
#     14
#     )

os.makedirs(output, exist_ok=True)

########################################################################
# Whether existing cached control
########################################################################

CACHEDIR = os.path.join(tempfile.gettempdir(), "DAJIN")
os.makedirs(CACHEDIR, exist_ok=True)

if cache_control.exists(control, CACHEDIR):
    IS_CACHED = True
else:
    cache_control.save_header(control, CACHEDIR)
    IS_CACHED = False

########################################################################
# minimap2
########################################################################
# importlib.reload(mapping)

mapping.split_fasta(allele, TMPDIR_PATHS["fasta"])

if IS_CACHED:
    cache_control.load(CACHEDIR, TMPDIR_PATHS["sam"])
else:
    mapping.minimap2(control, TMPDIR_PATHS["fasta"], TMPDIR_PATHS["sam"], threads)
    cache_control.save(control, TMPDIR_PATHS["sam"], CACHEDIR)

mapping.minimap2(sample, TMPDIR_PATHS["fasta"], TMPDIR_PATHS["sam"], threads)

########################################################################
# MIDS conversion
########################################################################
# importlib.reload(midsconv)

sampath = ".tmpDAJIN/sam/barcode31_albino.sam"
# SAMDIR = os.path.join("tests", "samTomids", "input")
# SAMFILES = [os.path.join(SAMDIR, _) for _ in os.listdir(SAMDIR)]

for samfile in os.listdir(TMPDIR_PATHS["sam"]):
    sampath = os.path.join(TMPDIR_PATHS["sam"], samfile)
    output = sampath.replace(".sam", ".csv").replace("sam", "midsconv")
    midscsv = midsconv.sam_to_mids(sampath, threads)
    with open(output, "w") as f:
        f.write("\n".join(midscsv))

########################################################################
# MIDS QC filtering
## 完全長リードのみを取り出す：両端から50bp連続して"="であるリードを除く
## Phread scoreが0.1以下のリードについて再分配する
########################################################################

print(mids_csv[0])

if debug is False:
    shutil.rmtree(".tmpDAJIN")

#! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

if __name__ == "__main__":
    main()
