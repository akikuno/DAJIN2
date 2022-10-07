from __future__ import annotations

import warnings

warnings.simplefilter("ignore")
import hashlib
from collections import defaultdict
from copy import deepcopy
from itertools import groupby
from pathlib import Path
from importlib import reload
import midsv
import pysam

from src.DAJIN2.core import preprocess, classification, clustering, consensus, report

reload(preprocess)
reload(classification)
reload(clustering)
reload(consensus)
reload(report)

# # * Point mutation
# SAMPLE, CONTROL, ALLELE, OUTPUT, GENOME, DEBUG, THREADS = (
#     "examples/pm-tyr/barcode31.fq.gz",
#     "examples/pm-tyr/barcode32.fq.gz",
#     "examples/pm-tyr/design_tyr.fa",
#     "DAJIN_results",
#     "mm10",
#     True,
#     14,
# )


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

# * flox insertion
SAMPLE, CONTROL, ALLELE, NAME, GENOME, DEBUG, THREADS = (
    "examples/flox-cables2/AyabeTask1/barcode36.fq.gz",
    "examples/flox-cables2/AyabeTask1/barcode42.fq.gz",
    "examples/flox-cables2/AyabeTask1/design_cables2.fa",
    "Ayabe-Task1",
    "mm10",
    True,
    14,
)

##########################################################
# Check inputs
##########################################################
preprocess.check_inputs.check_files(SAMPLE, CONTROL, ALLELE)
TEMPDIR = Path("DAJINResults", f".tempdir_{NAME}")
IS_CACHE_CONTROL = preprocess.check_inputs.is_cache_control(CONTROL, TEMPDIR)
IS_CACHE_GENOME = preprocess.check_inputs.is_cache_genome(GENOME, TEMPDIR, IS_CACHE_CONTROL)
UCSC_URL, GOLDENPATH_URL = None, None
if GENOME and not IS_CACHE_GENOME:
    UCSC_URL, GOLDENPATH_URL = preprocess.check_inputs.check_and_fetch_genome(GENOME)

##########################################################
# Format inputs
##########################################################
SAMPLE_NAME = preprocess.format_inputs.extract_basename(SAMPLE)
CONTROL_NAME = preprocess.format_inputs.extract_basename(CONTROL)
DICT_ALLELE = preprocess.format_inputs.dictionize_allele(ALLELE)

preprocess.format_inputs.make_directories(TEMPDIR)

if GENOME:
    GENOME_COODINATES, CHROME_SIZE = preprocess.format_inputs.get_coodinates_and_chromsize(
        TEMPDIR, GENOME, DICT_ALLELE, UCSC_URL, GOLDENPATH_URL, IS_CACHE_GENOME
    )


flag1 = Path(TEMPDIR, "midsv", f"{SAMPLE_NAME}_control.jsonl").exists()
flag2 = Path(TEMPDIR, "midsv", f"{CONTROL_NAME}_control.jsonl").exists()
flag = flag1 and flag2

if not flag:
    ################################################################################
    # Export fasta files as single-FASTA format
    ################################################################################
    # TODO: use yeild, not export
    for identifier, sequence in DICT_ALLELE.items():
        contents = "\n".join([">" + identifier, sequence]) + "\n"
        output_fasta = Path(TEMPDIR, "fasta", f"{identifier}.fasta")
        output_fasta.write_text(contents)
    ###############################################################################
    # Mapping with minimap2/mappy
    ###############################################################################
    for path_fasta in Path(TEMPDIR, "fasta").glob("*.fasta"):
        name_fasta = path_fasta.stem
        if name_fasta not in set(DICT_ALLELE.keys()):
            continue
        if not IS_CACHE_CONTROL:
            preprocess.mappy_align.output_sam(TEMPDIR, path_fasta, name_fasta, CONTROL, CONTROL_NAME)
        preprocess.mappy_align.output_sam(TEMPDIR, path_fasta, name_fasta, SAMPLE, SAMPLE_NAME)
    ########################################################################
    # MIDSV conversion
    ########################################################################
    if not IS_CACHE_CONTROL:
        for path_sam in Path(TEMPDIR, "sam").glob(f"{CONTROL_NAME}*"):
            preprocess.calc_midsv.output_midsv(TEMPDIR, path_sam, DICT_ALLELE)
    for path_sam in Path(TEMPDIR, "sam").glob(f"{SAMPLE_NAME}*"):
        preprocess.calc_midsv.output_midsv(TEMPDIR, path_sam, DICT_ALLELE)
    ###############################################################################
    # Cashe inputs (control)
    ###############################################################################
    if not IS_CACHE_CONTROL:
        control_hash = Path(CONTROL).read_bytes()
        control_hash = hashlib.sha256(control_hash).hexdigest()
        PATH_CACHE_HASH = Path(TEMPDIR, "cache", "control_hash.txt")
        PATH_CACHE_HASH.write_text(str(control_hash))

########################################################################
# Classify alleles
########################################################################
path_midsv = Path(TEMPDIR, "midsv").glob(f"{SAMPLE_NAME}*")
classif_sample = classification.classify_alleles(path_midsv, SAMPLE_NAME)

########################################################################
# Detect Structural variants
########################################################################
for classif in classif_sample:
    classif["SV"] = classification.detect_sv(classif["CSSPLIT"], threshold=50)

########################################################################
# Clustering
########################################################################
# -----------------------------------------------------------------------
# Extract significantly different base loci between Sample and Control
# -----------------------------------------------------------------------
dict_cssplit_control = defaultdict(list[dict])
for ALLELE in DICT_ALLELE.keys():
    path_control = Path(TEMPDIR, "midsv", f"{CONTROL_NAME}_{ALLELE}.jsonl")
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


# ALLELE = "control"
# SV = False
# cssplit_sample = []
# for group in classif_sample:
#     if group["ALLELE"] == ALLELE and group["SV"] == SV:
#         cssplit_sample.append(group["CSSPLIT"])

# table_sample, table_control = make_table(cssplit_sample, cssplit_control)
# s = table_sample[88]
# c = table_control[88]
# s, c
# -----------------------------------------------------------------------
# Clustering
# -----------------------------------------------------------------------
labels = []
label_start = 1
for (ALLELE, SV), group in groupby(classif_sample, key=lambda x: (x["ALLELE"], x["SV"])):
    key = f'{{"ALLELE": "{ALLELE}", "SV": {SV}}}'
    cssplit_sample = [g["CSSPLIT"] for g in group]
    diffloci = diffloci_by_alleles[key]
    if diffloci is None:
        scores = []
    else:
        scores = clustering.make_scores(cssplit_sample, diffloci)
    if any(scores):
        labels += [label + label_start for label in clustering.clustering(scores).tolist()]
    else:
        labels += [label_start] * len(cssplit_sample)
    label_start = len(set(labels)) + 1

clust_sample = deepcopy(classif_sample)
for clust, label in zip(clust_sample, labels):
    clust["LABEL"] = label
    del clust["CSSPLIT"]

n_sample = len(clust_sample)
d = defaultdict(int)
for cs in clust_sample:
    d[cs["LABEL"]] += 1 / n_sample

d_per = {key: round(val * 100, 1) for key, val in d.items()}

for cs in clust_sample:
    cs["PERCENT"] = d_per[cs["LABEL"]]

# Allocate new labels by PERCENT
clust_sample.sort(key=lambda x: (-x["PERCENT"], x["LABEL"]))
new_label = 1
prev_label = clust_sample[0]["LABEL"]
for cs in clust_sample:
    if prev_label != cs["LABEL"]:
        new_label += 1
    prev_label = cs["LABEL"]
    cs["LABEL"] = new_label

########################################################################
# Consensus call
########################################################################
path = Path(TEMPDIR, "midsv", f"{CONTROL_NAME}_control.jsonl")
cssplit_control = midsv.read_jsonl(path)
path = Path(TEMPDIR, "midsv", f"{SAMPLE_NAME}_control.jsonl")
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
    cons_fasta = report.report_files.to_fasta(header, cons_seq)
    Path(TEMPDIR, "reports", f"{SAMPLE_NAME}_{header}.fasta").write_text(cons_fasta)

# HTML
for header, cons_per in cons_percentage.items():
    cons_html = report.report_files.to_html(header, cons_per)
    Path(TEMPDIR, "reports", f"{SAMPLE_NAME}_{header}.html").write_text(cons_html)

# BAM
report.report_bam.output_bam(TEMPDIR, RESULT_SAMPLE, SAMPLE_NAME, GENOME, GENOME_COODINATES, CHROME_SIZE, THREADS)

# VCF
# working in progress


########################################################################
# Output Results
########################################################################
RESULT_SAMPLE = deepcopy(cssplit_sample)
for res in RESULT_SAMPLE:
    del res["RNAME"]
    del res["CSSPLIT"]

RESULT_SAMPLE.sort(key=lambda x: x["LABEL"])

d = defaultdict(int)
for res in RESULT_SAMPLE:
    d[res["NAME"]] += 1

d

for c in cons_percentage["allele6_flox_mutated"]:
    key = list(c.keys())[0]
    if "+" in key:
        print(list(c.items())[0])

