import shutil
import sys
import os
import tempfile

# Custom modules
import src.DAJIN2.utils.util as util
from src.DAJIN2.utils import cache_control
from src.DAJIN2.utils import argparser
from src.DAJIN2.preprocess import mids
from src.DAJIN2.preprocess import mapping

# For development
import importlib


#! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# def main():
########################################################################
# Setting
########################################################################

util.check_dependencies(["minimap2", "samtools"])

TMPDIR = ".tmpDAJIN"
os.makedirs(TMPDIR, exist_ok=True)
util.make_directories(TMPDIR, ["fasta", "sam"])
TMPDIR_PATHS = {_: os.path.join(TMPDIR, _) for _ in os.listdir(TMPDIR)}

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
# importlib.reload(mids)

sampath = ".tmpDAJIN/sam/barcode31_albino.sam"
# SAMDIR = os.path.join("tests", "samTomids", "input")
# SAMFILES = [os.path.join(SAMDIR, _) for _ in os.listdir(SAMDIR)]

mids_csv = mids.sam_to_mids(sampath)
print(mids_csv[0])

if debug is False:
    shutil.rmtree(".tmpDAJIN")

#! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

if __name__ == "__main__":
    main()
