from __future__ import annotations

import warnings

warnings.simplefilter("ignore")
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
    "examples/flox-cables2/AyabeTask1/barcode33.fq.gz",
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
TEMPDIR = Path("DAJINResults", ".tempdir", NAME)
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

MASKS_CONTROL = clustering.mask_control(TEMPDIR, DICT_ALLELE, CONTROL_NAME)

diffloci_by_alleles = clustering.extract_different_loci(
    TEMPDIR, classif_sample, MASKS_CONTROL, DICT_ALLELE, CONTROL_NAME
)
# ALLELE = "flox"
# SV = False
# cssplit_sample = [cs["CSSPLIT"] for cs in classif_sample if cs["ALLELE"] == ALLELE and cs["SV"] == SV]

clust_sample = clustering.add_labels(classif_sample, diffloci_by_alleles)

clust_sample = clustering.add_readnum(clust_sample)

clust_sample = clustering.add_percent(clust_sample)

clust_sample = clustering.update_labels(clust_sample)

########################################################################
# Consensus call
########################################################################

RESULT_SAMPLE, cons_percentage, cons_sequence = consensus.call(clust_sample, diffloci_by_alleles, DICT_ALLELE)


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

# BAM
report.report_bam.output_bam(TEMPDIR, RESULT_SAMPLE, SAMPLE_NAME, GENOME, GENOME_COODINATES, CHROME_SIZE, THREADS)

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

