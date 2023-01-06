from __future__ import annotations

import hashlib
from collections import defaultdict
from pathlib import Path
from importlib import reload

from src.DAJIN2.core import preprocess, classification, clustering, consensus, report

reload(preprocess)
reload(classification)
reload(clustering)
reload(consensus)
reload(report)

from src.DAJIN2.core.clustering import clustering

# # * Subset of Point mutation
# # 50 or 10 or 01%
# SAMPLE, CONTROL, ALLELE, NAME, GENOME, DEBUG, THREADS = (
#     "misc/data/tyr_albino_01%.fq.gz",
#     "misc/data/tyr_control.fq.gz",
#     "misc/data/tyr_control.fasta",
#     "test_tyr_albino_01%",
#     "mm10",
#     True,
#     14,
# )

# # * Point mutation
# SAMPLE, CONTROL, ALLELE, NAME, GENOME, DEBUG, THREADS = (
#     "examples/pm-tyr/barcode31.fq.gz",
#     "examples/pm-tyr/barcode32.fq.gz",
#     "examples/pm-tyr/design_tyr.fa",
#     "test-pm-tyr",
#     "mm10",
#     True,
#     14,
# )


# * 2-cut deletion
SAMPLE, CONTROL, ALLELE, NAME, GENOME, DEBUG, THREADS = (
    "tests/data/knockout/test_barcode25.fq.gz",
    "tests/data/knockout/test_barcode30.fq.gz",
    "tests/data/knockout/design_stx2.fa",
    "test-knockout",
    "mm10",
    True,
    14,
)

# #* 2-cut deletion
# SAMPLE, CONTROL, ALLELE, NAME, GENOME, DEBUG, THREADS = (
#     "examples/del-stx2/barcode25.fq.gz",
#     "examples/del-stx2/barcode30.fq.gz",
#     "examples/del-stx2/design_stx2.fa",
#     "test-stx2-deletion",
#     "mm10",
#     True,
#     14,
# )

# # * flox insertion
# SAMPLE, CONTROL, ALLELE, NAME, GENOME, DEBUG, THREADS = (
#     "examples/flox-cables2/AyabeTask1/barcode31.fq.gz",
#     "examples/flox-cables2/AyabeTask1/barcode42.fq.gz",
#     "examples/flox-cables2/AyabeTask1/design_cables2.fa",
#     "test-Ayabe-Task1",
#     "mm10",
#     True,
#     14,
# )

##########################################################
# Check inputs
##########################################################
preprocess.check_inputs.check_files(SAMPLE, CONTROL, ALLELE)
TEMPDIR = Path("DAJINResults", ".tempdir", NAME)
IS_CACHE_CONTROL = preprocess.check_inputs.exists_cached_control(CONTROL, TEMPDIR)
IS_CACHE_GENOME = preprocess.check_inputs.exists_cached_genome(GENOME, TEMPDIR, IS_CACHE_CONTROL)
UCSC_URL, GOLDENPATH_URL = None, None
if GENOME and not IS_CACHE_GENOME:
    UCSC_URL, GOLDENPATH_URL = preprocess.check_inputs.check_and_fetch_genome(GENOME)

##########################################################
# Format inputs
##########################################################
SAMPLE_NAME = preprocess.format_inputs.extract_basename(SAMPLE)
CONTROL_NAME = preprocess.format_inputs.extract_basename(CONTROL)
FASTA_ALLELES = preprocess.format_inputs.dictionize_allele(ALLELE)

preprocess.format_inputs.make_directories(TEMPDIR, SAMPLE_NAME, CONTROL_NAME)

if GENOME:
    GENOME_COODINATES = preprocess.format_inputs.fetch_coodinate(GENOME, UCSC_URL, FASTA_ALLELES["control"])
    CHROME_SIZE = preprocess.format_inputs.fetch_chrom_size(GENOME_COODINATES["chr"], GENOME, GOLDENPATH_URL)
    preprocess.format_inputs.cache_coodinates_and_chromsize(TEMPDIR, GENOME, GENOME_COODINATES, CHROME_SIZE)

flag1 = Path(TEMPDIR, "midsv", f"{CONTROL_NAME}_splice_control.jsonl").exists()
flag2 = Path(TEMPDIR, "midsv", f"{SAMPLE_NAME}_splice_control.jsonl").exists()
flag = flag1 and flag2

if not flag:
    ################################################################################
    # Export fasta files as single-FASTA format
    ################################################################################
    # TODO: use yeild, not export
    for identifier, sequence in FASTA_ALLELES.items():
        contents = "\n".join([">" + identifier, sequence]) + "\n"
        output_fasta = Path(TEMPDIR, "fasta", f"{identifier}.fasta")
        output_fasta.write_text(contents)
    ###############################################################################
    # Mapping with minimap2/mappy
    ###############################################################################
    for path_fasta in Path(TEMPDIR, "fasta").glob("*.fasta"):
        name_fasta = path_fasta.stem
        preprocess.mappy_align.output_sam(TEMPDIR, path_fasta, name_fasta, CONTROL, CONTROL_NAME, threads=THREADS)
        preprocess.mappy_align.output_sam(TEMPDIR, path_fasta, name_fasta, SAMPLE, SAMPLE_NAME, threads=THREADS)
        preprocess.mappy_align.output_sam(
            TEMPDIR, path_fasta, name_fasta, CONTROL, CONTROL_NAME, preset="splice", threads=THREADS
        )
        preprocess.mappy_align.output_sam(
            TEMPDIR, path_fasta, name_fasta, SAMPLE, SAMPLE_NAME, preset="splice", threads=THREADS
        )
    ########################################################################
    # MIDSV conversion
    ########################################################################
    for path_sam in Path(TEMPDIR, "sam").glob(f"{CONTROL_NAME}_splice_*"):
        preprocess.calc_midsv.output_midsv(TEMPDIR, path_sam)
    for path_sam in Path(TEMPDIR, "sam").glob(f"{SAMPLE_NAME}_splice_*"):
        preprocess.calc_midsv.output_midsv(TEMPDIR, path_sam)
    ###############################################################################
    # Correct CSSPLITS
    ###############################################################################
    preprocess.correct_revititive_deletions.execute(TEMPDIR, FASTA_ALLELES, CONTROL_NAME, SAMPLE_NAME)
    preprocess.correct_sequence_error.execute(TEMPDIR, FASTA_ALLELES, CONTROL_NAME, SAMPLE_NAME)
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

paths_midsv = list(Path(TEMPDIR, "midsv").glob(f"{SAMPLE_NAME}_splice*"))
classif_sample = classification.classify_alleles(paths_midsv)

# paths_midsv = list(Path(TEMPDIR, "midsv").glob(f"{CONTROL_NAME}*"))
# classif_control = classification.classify_alleles(paths_midsv)

########################################################################
# Detect Structural variants
########################################################################

for classif in classif_sample:
    classif["SV"] = classification.detect_sv(classif["CSSPLIT"], threshold=50)

# for classif in classif_control:
#     classif["SV"] = classification.detect_sv(classif["CSSPLIT"], threshold=50)

# d = defaultdict(int)
# for c in classif_sample:
#     keys = f'{c["ALLELE"]}-{c["SV"]}-{c["PRESET"]}'
#     d[keys] += 1

# d

########################################################################
# Clustering
########################################################################

# allele = "control"
# sv = False
# cssplit_control = [cs["CSSPLIT"] for cs in classif_control if cs["ALLELE"] == allele and cs["SV"] == sv]

# for classif in classif_sample:
#     if "e0174c91bd7b" in classif["QNAME"]:
#         print(classif["CSSPLIT"])

# # 476 del

# for classif in classif_sample:
#     if "edc66c43a83f" in classif["QNAME"]:
#         print(classif["CSSPLIT"])

# 615 del

# KNOCKIN_LOCI = clustering.find_knockin_loci(TEMPDIR, FASTA_ALLELES, CONTROL_NAME)

# DIFFLOCI_ALLELES, REPETITIVE_DELLOCI = clustering.extract_different_loci(
#     TEMPDIR, classif_sample, KNOCKIN_LOCI, FASTA_ALLELES, CONTROL_NAME
# )

clust_sample = clustering.add_labels(classif_sample, TEMPDIR, CONTROL_NAME, FASTA_ALLELES, THREADS)
clust_sample = clustering.add_readnum(clust_sample)
clust_sample = clustering.add_percent(clust_sample)
clust_sample = clustering.update_labels(clust_sample)

# import midsv
# from pathlib import Path
# Save
# # midsv.write_jsonl(clust_sample, Path(TEMPDIR, "clustering", f"{SAMPLE_NAME}.jsonl"))
# # Load
# # clust_sample = midsv.read_jsonl(Path(TEMPDIR, "clustering", f"{SAMPLE_NAME}.jsonl"))

# d_count = defaultdict(int)
# for c in clust_sample:
#     d_count[c["LABEL"]] += 1

# coverage = len(clust_sample)
# d_percent = defaultdict(int)
# for c in clust_sample:
#     d_percent[c["LABEL"]] += 1 / coverage * 100

# d_count
# d_percent

# d = defaultdict(int)
# for c in clust_sample:
#     keys = f'{c["ALLELE"]}-{c["SV"]}-{c["LABEL"]}'
#     d[keys] += 1

# d
# from pprint import pprint

# pprint(d)
########################################################################
# Consensus call
########################################################################

RESULT_SAMPLE, cons_percentage, cons_sequence = consensus.call(clust_sample, FASTA_ALLELES)
# cons_percentage["allele1_flox_mutated_52.968%"][461]
# keys = []
# keys.append("flox")
# keys.append(False)
# cssplit_sample = [cs for cs in clust_sample if cs["ALLELE"] == keys[0] and cs["SV"] == keys[1]]

# readid = "3f7aff2243b8"  # strand +
# readid = "28e12f325a68"  # strand -
# for res in RESULT_SAMPLE:
#     if readid in res["QNAME"]:
#         print(res)

# ----------------------------------------------------------
# Conseusns Report：FASTA/HTML/VCF
# ----------------------------------------------------------
# FASTA
for header, cons_seq in cons_sequence.items():
    cons_fasta = report.report_files.to_fasta(header, cons_seq)
    Path(TEMPDIR, "report", "FASTA", f"{SAMPLE_NAME}_{header}.fasta").write_text(cons_fasta)

# HTML
for header, cons_per in cons_percentage.items():
    cons_html = report.report_files.to_html(header, cons_per)
    Path(TEMPDIR, "report", "HTML", f"{SAMPLE_NAME}_{header}.html").write_text(cons_html)

# BAM and igvjs
report.report_bam.output_bam_control(TEMPDIR, CONTROL_NAME, GENOME, GENOME_COODINATES, CHROME_SIZE, THREADS)
report.report_bam.output_bam_sample(
    TEMPDIR, RESULT_SAMPLE, SAMPLE_NAME, GENOME, GENOME_COODINATES, CHROME_SIZE, THREADS
)

# VCF
# working in progress


########################################################################
# MEMO
########################################################################

d = defaultdict(int)
for res in RESULT_SAMPLE:
    d[res["NAME"]] += 1

d

# ALLELE = "control"
# SV = False
# cssplit_sample = []
# for group in classif_sample:
#     if group["ALLELE"] == ALLELE and group["SV"] == SV:
#         cssplit_sample.append(group["CSSPLIT"])
