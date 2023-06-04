from __future__ import annotations

import sys
import os
from pathlib import Path

# 初期は必ずDAJIN2-experiments/以下の階層にあるので、DAJIN2-experimentsまで戻る
path_parent = Path.cwd()
while path_parent.stem != "DAJIN2-experiments":
    path_parent = path_parent.parent

# DAJIN2-experiments/DAJIN2のディレクトリに移動し、PATHを通す
PATH_DAJIN2 = Path.joinpath(path_parent, "DAJIN2")
os.chdir(PATH_DAJIN2)
sys.path.append(PATH_DAJIN2)
sys.path.append(Path.joinpath(PATH_DAJIN2, "src"))


import hashlib
from collections import defaultdict
from importlib import reload

from src.DAJIN2.core import preprocess, classification, clustering, consensus, report

reload(preprocess)
reload(classification)
reload(clustering)
reload(consensus)
reload(report)

# #### # * Subset of Point mutation
# #### # 50 or 10 or 01%
# percent = "01"
# SAMPLE, CONTROL, ALLELE, NAME, GENOME, DEBUG, THREADS = (
#     f"misc/data/tyr_albino_{percent}%.fq.gz",
#     "misc/data/tyr_control.fq.gz",
#     "misc/data/tyr_control.fasta",
#     f"test-tyr-albino-{percent}%",
#     "mm10",
#     True,
#     30,
# )

# # ##### # * R9 vs R10 Point mutation
# SAMPLE, CONTROL, ALLELE, NAME, GENOME, DEBUG, THREADS = (
#     "../fastq_20230304/barcode31.fq.gz",
#     "../fastq_20230304/barcode32.fq.gz",
#     "examples/pm-tyr/design_tyr.fa",
#     "test-20230304-pm-tyr",
#     "mm10",
#     True,
#     6,
# )

# ##### # * Point mutation
# SAMPLE, CONTROL, ALLELE, NAME, GENOME, DEBUG, THREADS = (
#     "examples/pm-tyr/barcode31.fq.gz",
#     "examples/pm-tyr/barcode32.fq.gz",
#     "examples/pm-tyr/design_tyr.fa",
#     "test-pm-tyr",
#     "mm10",
#     True,
#     14,
# )


# ##### # * 2-cut deletion
# SAMPLE, CONTROL, ALLELE, NAME, GENOME, DEBUG, THREADS = (
#     "tests/data/knockout/test_barcode25.fq.gz",
#     "tests/data/knockout/test_barcode30.fq.gz",
#     "tests/data/knockout/design_stx2.fa",
#     "single-stx2deletion",
#     "mm10",
#     True,
#     30,
# )

# #### #* 2-cut deletion
SAMPLE, CONTROL, ALLELE, NAME, GENOME, DEBUG, THREADS = (
    "examples/del-stx2/barcode25.fq.gz",
    "examples/del-stx2/barcode30.fq.gz",
    "examples/del-stx2/design_stx2.fa",
    "test-stx2-deletion",
    "mm10",
    True,
    14,
)

# #### * flox insertion
# SAMPLE, CONTROL, ALLELE, NAME, GENOME, DEBUG, THREADS = (
#     "examples/flox-cables2/AyabeTask1/barcode31.fq.gz",
#     "examples/flox-cables2/AyabeTask1/barcode42.fq.gz",
#     "examples/flox-cables2/AyabeTask1/design_cables2.fa",
#     "single-ayabetask1",
#     "mm10",
#     True,
#     14,
# )

######################################################################
# Preprocessing
######################################################################

print(f"processing {NAME}...")

SAMPLE = preprocess.format_inputs.convert_to_posix_path(SAMPLE)
CONTROL = preprocess.format_inputs.convert_to_posix_path(CONTROL)
ALLELE = preprocess.format_inputs.convert_to_posix_path(ALLELE)

# ====================================================================
# Varidate inputs
# ====================================================================

preprocess.validate_inputs.check_files(SAMPLE, CONTROL, ALLELE)
TEMPDIR = Path("DAJINResults", ".tempdir", NAME)
IS_CACHE_CONTROL = preprocess.validate_inputs.exists_cached_control(CONTROL, TEMPDIR)
IS_CACHE_GENOME = preprocess.validate_inputs.exists_cached_genome(GENOME, TEMPDIR, IS_CACHE_CONTROL)
UCSC_URL, GOLDENPATH_URL = None, None
if GENOME and not IS_CACHE_GENOME:
    UCSC_URL, GOLDENPATH_URL = preprocess.validate_inputs.check_and_fetch_genome(GENOME)

# ====================================================================
# Format inputs
# ====================================================================
SAMPLE_NAME = preprocess.format_inputs.extract_basename(SAMPLE)
CONTROL_NAME = preprocess.format_inputs.extract_basename(CONTROL)
FASTA_ALLELES = preprocess.format_inputs.dictionize_allele(ALLELE)
THREADS = preprocess.format_inputs.update_threads(THREADS)

SUBDIRS = ["cache", "fasta", "sam", "midsv", "midsv_corrected", "report", "result", "mutation_loci"]
preprocess.format_inputs.make_directories(TEMPDIR, SUBDIRS, SAMPLE_NAME, CONTROL_NAME)

if GENOME:
    GENOME_COODINATES = preprocess.format_inputs.fetch_coodinate(GENOME, UCSC_URL, FASTA_ALLELES["control"])
    CHROME_SIZE = preprocess.format_inputs.fetch_chrom_size(GENOME_COODINATES["chr"], GENOME, GOLDENPATH_URL)
    preprocess.format_inputs.cache_coodinates_and_chromsize(TEMPDIR, GENOME, GENOME_COODINATES, CHROME_SIZE)

flag1 = Path(TEMPDIR, "midsv", f"{CONTROL_NAME}_splice_control.jsonl").exists()
flag2 = Path(TEMPDIR, "midsv", f"{SAMPLE_NAME}_splice_control.jsonl").exists()
flag = flag1 and flag2

if not flag:
    # ====================================================================
    # Export fasta files as single-FASTA format
    # ====================================================================
    # TODO: use yeild, not export
    for identifier, sequence in FASTA_ALLELES.items():
        contents = "\n".join([">" + identifier, sequence]) + "\n"
        output_fasta = Path(TEMPDIR, "fasta", f"{identifier}.fasta")
        output_fasta.write_text(contents)
    # ====================================================================
    # Mapping with mappy
    # ====================================================================
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
    # ====================================================================
    # MIDSV conversion
    # ====================================================================
    for allele in FASTA_ALLELES:
        preprocess.call_midsv(TEMPDIR, CONTROL_NAME, allele)
    for allele in FASTA_ALLELES:
        preprocess.call_midsv(TEMPDIR, SAMPLE_NAME, allele)
    # ====================================================================
    # Convert any `N` as deletions other than consecutive `N` from both ends
    # ====================================================================
    preprocess.replace_NtoD(TEMPDIR, FASTA_ALLELES, CONTROL_NAME)
    preprocess.replace_NtoD(TEMPDIR, FASTA_ALLELES, SAMPLE_NAME)
    # ====================================================================
    # CSSPLITS Error Correction
    # ====================================================================
    MUTATION_LOCI_ALLELES = preprocess.extract_mutation_loci(TEMPDIR, FASTA_ALLELES, CONTROL_NAME, SAMPLE_NAME)
    preprocess.correct_sequence_error(TEMPDIR, FASTA_ALLELES, CONTROL_NAME, SAMPLE_NAME, MUTATION_LOCI_ALLELES)
    # preprocess.correct_knockin.execute(TEMPDIR, FASTA_ALLELES, CONTROL_NAME, SAMPLE_NAME)
    # ====================================================================
    # Cashe inputs (control)
    # ====================================================================
    if not IS_CACHE_CONTROL:
        control_hash = Path(CONTROL).read_bytes()
        control_hash = hashlib.sha256(control_hash).hexdigest()
        PATH_CACHE_HASH = Path(TEMPDIR, "cache", "control_hash.txt")
        PATH_CACHE_HASH.write_text(str(control_hash))

####################################################################################
# Classify alleles
####################################################################################
print("Classify...")

classif_sample = classification.classify_alleles(TEMPDIR, SAMPLE_NAME)

# for classif in classif_sample:
#     classif["SV"] = classification.detect_sv(classif["CSSPLIT"], threshold=50)

####################################################################################
# Clustering
####################################################################################
print("Clustering...")

MUTATION_LOCI = clustering.extract_mutation_loci(TEMPDIR, FASTA_ALLELES, SAMPLE_NAME, CONTROL_NAME)
clust_sample = clustering.add_labels(classif_sample, TEMPDIR, CONTROL_NAME, MUTATION_LOCI, THREADS)
clust_sample = clustering.add_readnum(clust_sample)
clust_sample = clustering.add_percent(clust_sample)
clust_sample = clustering.update_labels(clust_sample)

####################################################################################
# Consensus call
####################################################################################
print("Consensus call...")

cons_percentage, cons_sequence = consensus.call_consensus(clust_sample)
allele_names = consensus.call_allele_name(cons_sequence, cons_percentage, FASTA_ALLELES)
cons_percentage = consensus.update_key_by_allele_name(cons_percentage, allele_names)
cons_sequence = consensus.update_key_by_allele_name(cons_sequence, allele_names)
RESULT_SAMPLE = consensus.add_key_by_allele_name(clust_sample, allele_names)
RESULT_SAMPLE.sort(key=lambda x: x["LABEL"])

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


####################################################################################
# MEMO
####################################################################################

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
