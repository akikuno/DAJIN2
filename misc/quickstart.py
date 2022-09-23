from __future__ import annotations

import warnings

warnings.simplefilter("ignore")
import shutil
from collections import defaultdict
from copy import deepcopy
from itertools import groupby
from pathlib import Path
from urllib.error import URLError

import midsv
import pysam

from src.DAJIN2 import classification, clustering
from src.DAJIN2.consensus import module_consensus as consensus
from src.DAJIN2.preprocess import argparser, check_inputs, format_inputs, mappy_align
from src.DAJIN2.report import report_af, report_bam

# * Point mutation
SAMPLE, CONTROL, ALLELE, OUTPUT, GENOME, DEBUG, THREADS = (
    "examples/pm-tyr/barcode31.fq.gz",
    "examples/pm-tyr/barcode32.fq.gz",
    "examples/pm-tyr/design_tyr.fa",
    "DAJIN_results",
    "mm10",
    True,
    14,
)


# # * 2-cut deletion
# SAMPLE, CONTROL, ALLELE, OUTPUT, GENOME, DEBUG, THREADS = (
#     "examples/del-stx2/barcode25.fq.gz",
#     "examples/del-stx2/barcode30.fq.gz",
#     "examples/del-stx2/design_stx2.fa",
#     "DAJIN_results",
#     "mm10",
#     True,
#     14,
# )

# # * flox insertion
# SAMPLE, CONTROL, ALLELE, OUTPUT, GENOME, DEBUG, THREADS = (
#     "examples/flox-cables2/AyabeTask1/barcode31.fq.gz",
#     "examples/flox-cables2/AyabeTask1/barcode42.fq.gz",
#     "examples/flox-cables2/AyabeTask1/design_cables2.fa",
#     "DAJIN_results",
#     "mm10",
#     True,
#     14,
# )


###############################################################################
# Check and format inputs (sample/control/allele)
###############################################################################

from pathlib import Path
from importlib import reload
from src.DAJIN2.preprocess import check_inputs
from src.DAJIN2.preprocess import format_inputs
from urllib.error import URLError

reload(check_inputs)
reload(format_inputs)

###############################################################################
# Check inputs (SAMPLE/CONTROL/ALLELE/GENOME)
###############################################################################

# ------------------------------------------------------------------------------
# Check input path
# ------------------------------------------------------------------------------

check_inputs.exists(CONTROL)
check_inputs.exists(SAMPLE)
check_inputs.exists(ALLELE)

# ------------------------------------------------------------------------------
# Check formats (extensions and contents)
# ------------------------------------------------------------------------------

check_inputs.fastq_extension(CONTROL)
check_inputs.fastq_content(CONTROL)
check_inputs.fastq_extension(SAMPLE)
check_inputs.fastq_content(SAMPLE)
check_inputs.fasta_content(ALLELE)

# ------------------------------------------------------------------------------
# Check GENOMEs if GENOME is inputted
# ------------------------------------------------------------------------------

if GENOME:
    # Check UCSC Server
    UCSC_URLS = [
        "https://genome.ucsc.edu/",
        "https://genome-asia.ucsc.edu/",
        "https://genome-euro.ucsc.edu/",
    ]
    UCSC_URL, flag_fail = check_inputs.available_url(UCSC_URLS)
    if flag_fail:
        raise URLError("UCSC Servers are currently down")
    # Check UCSC Download Server
    GOLDENPATH_URLS = [
        "https://hgdownload.cse.ucsc.edu/goldenPath",
        "http://hgdownload-euro.soe.ucsc.edu/goldenPath",
    ]
    GOLDENPATH_URL, flag_fail = check_inputs.available_url(GOLDENPATH_URLS)
    if flag_fail:
        raise URLError("UCSC Download Servers are currently down")
    # Check input genome
    check_inputs.available_genome(GENOME, UCSC_URL)


###############################################################################
# Format inputs (SAMPLE/CONTROL/ALLELE/GENOME)
###############################################################################

SAMPLE_NAME = format_inputs.extract_basename(SAMPLE)
CONTROL_NAME = format_inputs.extract_basename(CONTROL)
DICT_ALLELE = format_inputs.dictionize_allele(ALLELE)

if GENOME:
    GENOME_COODINATES = format_inputs.fetch_coodinate(GENOME, UCSC_URL, DICT_ALLELE["control"])
    CHROME_SIZE = format_inputs.fetch_chrom_size(GENOME_COODINATES["chr"], GENOME, GOLDENPATH_URL)

flag1 = Path(OUTPUT, ".tempdir", "midsv", f"{SAMPLE_NAME}_control.jsonl").exists()
flag2 = Path(OUTPUT, ".tempdir", "midsv", f"{CONTROL_NAME}_control.jsonl").exists()
flag = flag1 and flag2

if not flag:
    if Path(OUTPUT).exists():
        shutil.rmtree(OUTPUT)
    Path(OUTPUT).mkdir(exist_ok=True)
    if Path(OUTPUT, ".tempdir").exists():
        shutil.rmtree(Path(OUTPUT, ".tempdir"))
    for subdir in ["fasta", "fastq", "sam", "midsv", "bam", "reports"]:
        Path(OUTPUT, ".tempdir", subdir).mkdir(parents=True, exist_ok=True)
    # TODO: use yeild, not export
    for header, sequence in DICT_ALLELE.items():
        contents = "\n".join([">" + header, sequence]) + "\n"
        output_fasta = Path(OUTPUT, ".tempdir", "fasta", f"{header}.fasta")
        output_fasta.write_text(contents)
    ###############################################################################
    # Mapping with minimap2/mappy
    ###############################################################################
    for input_fasta in Path(OUTPUT, ".tempdir", "fasta").glob("*.fasta"):
        fasta_name = input_fasta.name.replace(".fasta", "")
        for fastq, fastq_name in zip([CONTROL, SAMPLE], [CONTROL_NAME, SAMPLE_NAME]):
            # Todo: 並行処理で高速化！！
            sam = mappy_align.to_sam(str(input_fasta), fastq)
            output_sam = Path(OUTPUT, ".tempdir", "sam", f"{fastq_name}_{fasta_name}.sam")
            output_sam.write_text("\n".join(sam))
    ########################################################################
    # MIDSV conversion
    ########################################################################
    for sampath in Path(OUTPUT, ".tempdir", "sam").iterdir():
        output_jsonl = Path(OUTPUT, ".tempdir", "midsv", f"{sampath.stem}.jsonl")
        sam = midsv.read_sam(sampath)
        midsv_jsonl = midsv.transform(sam, midsv=False, cssplit=True, qscore=False)
        midsv.write_jsonl(midsv_jsonl, output_jsonl)

########################################################################
# Classify alleles
########################################################################


path_midsv = Path(OUTPUT, ".tempdir", "midsv").glob(f"{SAMPLE_NAME}*")
classif_sample = classification.classify_alleles(path_midsv, SAMPLE_NAME)

path_midsv = Path(OUTPUT, ".tempdir", "midsv").glob(f"{CONTROL_NAME}*")
classif_control = classification.classify_alleles(path_midsv, CONTROL_NAME)

########################################################################
# Detect Structural variants
########################################################################

for classifs in [classif_sample, classif_control]:
    for classif in classifs:
        classif["SV"] = classification.detect_sv(classif["CSSPLIT"], threshold=50)

########################################################################
# Clustering
########################################################################

# -----------------------------------------------------------------------
# Extract significantly different base loci between Sample and Control
# -----------------------------------------------------------------------


dict_cssplit_control = defaultdict(list[dict])
for ALLELE in DICT_ALLELE.keys():
    path_control = Path(OUTPUT, ".tempdir", "midsv", f"{CONTROL_NAME}_{ALLELE}.jsonl")
    cssplit_control = [cs["CSSPLIT"] for cs in midsv.read_jsonl(path_control)]
    dict_cssplit_control[ALLELE] = cssplit_control

classif_sample.sort(key=lambda x: (x["ALLELE"], x["SV"]))
diffloci_by_alleles = defaultdict(list[dict])
for (ALLELE, SV), group in groupby(classif_sample, key=lambda x: (x["ALLELE"], x["SV"])):
    cssplit_sample = [record["CSSPLIT"] for record in group]
    cssplit_control = dict_cssplit_control[ALLELE]
    sequence = DICT_ALLELE[ALLELE]
    diffloci = clustering.screen_different_loci(cssplit_sample, cssplit_control, sequence, alpha=0.01, threshold=0.05)
    diffloci_by_alleles[f'{{"ALLELE": "{ALLELE}", "SV": {SV}}}'] = diffloci

# -----------------------------------------------------------------------
# Clustering
# -----------------------------------------------------------------------

labels = []
label_start = 1
for (ALLELE, SV), group in groupby(classif_sample, key=lambda x: (x["ALLELE"], x["SV"])):
    key = f'{{"ALLELE": "{ALLELE}", "SV": {SV}}}'
    cssplit_sample = [g["CSSPLIT"] for g in group]
    diffloci = diffloci_by_alleles[key]
    scores = list(clustering.make_scores(cssplit_sample, diffloci))
    if any(scores):
        labels += [label + label_start for label in clustering.clustering(scores).tolist()]
    else:
        labels += [label_start] * len(cssplit_sample)
    label_start = len(set(labels)) + 1

clust_sample = deepcopy(classif_sample)
for clust, label in zip(clust_sample, labels):
    clust["LABEL"] = label
    del clust["CSSPLIT"]

clust_sample.sort(key=lambda x: (x["ALLELE"], x["SV"], x["LABEL"]))

########################################################################
# Consensus call
########################################################################

path = Path(OUTPUT, ".tempdir", "midsv", f"{CONTROL_NAME}_control.jsonl")
cssplit_control = midsv.read_jsonl(path)

path = Path(OUTPUT, ".tempdir", "midsv", f"{SAMPLE_NAME}_control.jsonl")
cssplit_sample = midsv.read_jsonl(path)
cssplit_sample = consensus.join_listdicts(clust_sample, cssplit_sample, key="QNAME")
cssplit_sample.sort(key=lambda x: (x["ALLELE"], x["SV"], x["LABEL"]))

cons_percentage = defaultdict(list)
cons_sequence = defaultdict(list)
for keys, cssplits in groupby(cssplit_sample, key=lambda x: (x["ALLELE"], x["SV"], x["LABEL"])):
    cssplits = list(cssplits)
    cons_per = consensus.call_percentage(cssplits, cssplit_control)
    cons_seq = consensus.call_sequence(cons_per)
    allele_name = consensus.call_allele_name(keys, cons_seq, DICT_ALLELE)
    cons_percentage[allele_name] = cons_per
    cons_sequence[allele_name] = cons_seq
    for cs in cssplits:
        cs["NAME"] = allele_name

result_sample = deepcopy(cssplit_sample)
for res in result_sample:
    del res["RNAME"]
    del res["CSSPLIT"]

