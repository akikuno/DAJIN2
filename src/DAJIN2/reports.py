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

#! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# def main():

########################################################################
# Argument parse
########################################################################

SAMPLE, CONTROL, ALLELE, OUTPUT, GENOME, DEBUG, THREADS = argparser.parse()

########################################################################
# Whether existing cached CONTROL
########################################################################

# CACHEDIR = os.path.join(tempfile.gettempdir(), "DAJIN")
# os.makedirs(CACHEDIR, exist_ok=True)

# if cache_CONTROL.exists(CONTROL, CACHEDIR):
#     IS_CACHED = True
# else:
#     cache_CONTROL.save_header(CONTROL, CACHEDIR)
#     IS_CACHED = False

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


########################################################################
# Make directories
########################################################################

if Path(OUTPUT).exists():
    shutil.rmtree(OUTPUT)

Path(OUTPUT).mkdir(exist_ok=True)

if Path(OUTPUT, ".tempdir").exists():
    shutil.rmtree(Path(OUTPUT, ".tempdir"))

subdirectoris = ["fasta", "sam", "midsv", "bam", "reports"]
for subdir in subdirectoris:
    Path(OUTPUT, ".tempdir", subdir).mkdir(parents=True, exist_ok=True)

###############################################################################
# Format inputs (SAMPLE/CONTROL/ALLELE/GENOME)
###############################################################################

SAMPLE_NAME = format_inputs.extract_basename(SAMPLE)
CONTROL_NAME = format_inputs.extract_basename(CONTROL)
DICT_ALLELE = format_inputs.dictionize_allele(ALLELE)

if GENOME:
    GENOME_COODINATES = format_inputs.fetch_coodinate(GENOME, UCSC_URL, DICT_ALLELE["control"])
    CHROME_SIZE = format_inputs.fetch_chrom_size(GENOME_COODINATES["chr"], GENOME, GOLDENPATH_URL)

################################################################################
# Export fasta files as single-FASTA format
################################################################################

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
        # Todo: Speed up by parallele processing
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

RESULT_SAMPLE = deepcopy(cssplit_sample)
for res in RESULT_SAMPLE:
    del res["RNAME"]
    del res["CSSPLIT"]


# ----------------------------------------------------------
# Conseusns Report：FASTA/HTML/VCF
# ----------------------------------------------------------

# FASTA
for header, cons_seq in cons_sequence.items():
    cons_fasta = consensus.to_fasta(header, cons_seq)
    Path(f"{OUTPUT}/.tempdir/reports/{SAMPLE_NAME}_{header}.fasta").write_text(cons_fasta)

# HTML
for header, cons_per in cons_percentage.items():
    cons_html = consensus.to_html(header, cons_per)
    Path(f"{OUTPUT}/.tempdir/reports/{SAMPLE_NAME}_{header}.html").write_text(cons_html)


# VCF
# working in progress

########################################################################
# Report：アレル割合
# sample, allele name, #read, %read
########################################################################

# ----------------------------------------------------------
# All data
# ----------------------------------------------------------

df_clust_sample = report_af.all_allele(RESULT_SAMPLE, SAMPLE_NAME)
df_clust_sample.to_csv(f"{OUTPUT}/.tempdir/reports/{SAMPLE_NAME}_read_classification_all.csv", index=False)
df_clust_sample.to_excel(f"{OUTPUT}/.tempdir/reports/{SAMPLE_NAME}_read_classification_all.xlsx", index=False)

# ----------------------------------------------------------
# Summary data
# ----------------------------------------------------------

df_allele_frequency = report_af.summary_allele(RESULT_SAMPLE, SAMPLE_NAME)
df_allele_frequency.to_csv(f"{OUTPUT}/.tempdir/reports/{SAMPLE_NAME}_read_classification_summary.csv", index=False)
df_allele_frequency.to_excel(f"{OUTPUT}/.tempdir/reports/{SAMPLE_NAME}_read_classification_summary.xlsx", index=False)

# ----------------------------------------------------------
# Visualization
# ----------------------------------------------------------
g = report_af.plot(df_allele_frequency)
g.save(filename=f"{OUTPUT}/.tempdir/reports/tmp_output.png", dpi=350)
g.save(filename=f"{OUTPUT}/.tempdir/reports/tmp_output.pdf")


########################################################################
# Report：BAM
########################################################################


Path(OUTPUT, ".tempdir", "bam", "tmp_sam").mkdir(parents=True, exist_ok=True)

path_sam = Path(OUTPUT, ".tempdir", "sam", f"{SAMPLE_NAME}_control.sam")
sam = midsv.read_sam(path_sam)

sam_update = sam.copy()
sam_update = report_bam.remove_overlapped_reads(sam_update)
sam_update = report_bam.remove_microhomology(sam_update)

path_sam = f"{OUTPUT}/.tempdir/bam/tmp_sam/{SAMPLE_NAME}_control_updated.sam"
if GENOME:
    sam_update = report_bam.realign(sam_update, GENOME_COODINATES, CHROME_SIZE)
    report_bam.write_sam(sam_update, path_sam)
else:
    report_bam.write_sam(sam_update, path_sam)

pysam.sort("-@", f"{THREADS}", "-o", f"{OUTPUT}/.tempdir/bam/{SAMPLE_NAME}.bam", path_sam)
pysam.index("-@", f"{THREADS}", f"{OUTPUT}/.tempdir/bam/{SAMPLE_NAME}.bam")

sam_headers = [s for s in sam_update if s[0].startswith("@")]
sam_contents = [s for s in sam_update if not s[0].startswith("@")]
sam_groups = report_bam.group_by_name(sam_contents, clust_sample)
for key, sam_content in sam_groups.items():
    path_sam = f"{OUTPUT}/.tempdir/bam/tmp_sam/{key}.sam"
    report_bam.write_sam(sam_headers + sam_content, path_sam)
    pysam.sort("-@", f"{THREADS}", "-o", f"{OUTPUT}/.tempdir/bam/{SAMPLE_NAME}_{key}.bam", path_sam)
    pysam.index("-@", f"{THREADS}", f"{OUTPUT}/.tempdir/bam/{SAMPLE_NAME}_{key}.bam")


########################################################################
# Report：IGV.js
########################################################################

# IGV.js（10本のリードのみ表示）のために各サンプル50本程度のリードのみを抽出する

# for path_sam in Path(OUTPUT, ".tempdir", "sam").glob("*_control.sam"):
#     sam = midsv.read_sam(path_sam)
#     pysam.sort("-o", f"{OUTPUT}/.tempdir/bam/{name}.bam", str(path_sam))
#     pysam.index(f"{OUTPUT}/.tempdir/bam/{name}.bam")


########################################################################
# Finish call
########################################################################

if not DEBUG:
    shutil.rmtree(Path(OUTPUT, ".tempdir"))

# if __name__ == "__main__":
#     main()
