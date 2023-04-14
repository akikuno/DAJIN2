from __future__ import annotations

import hashlib
from collections import defaultdict
from pathlib import Path
from importlib import reload
import shutil
from src.DAJIN2.core import preprocess, classification, clustering, consensus, postprocess

reload(preprocess)
reload(classification)
reload(clustering)
reload(consensus)
reload(postprocess)

# * Point mutation
name = "tyr_albino_50%_only_control"
shutil.rmtree(f"DAJINResults/.tempdir/{name}/report/")
SAMPLE, CONTROL, ALLELE, NAME, GENOME, DEBUG, THREADS = (
    f"misc/data/tyr_albino_50%.fq.gz",
    "misc/data/tyr_control.fq.gz",
    "misc/data/tyr_control.fasta",
    f"{name}",
    "mm10",
    True,
    14,
)

##########################################################
# Check inputs
##########################################################
preprocess.validate_inputs.check_files(SAMPLE, CONTROL, ALLELE)
TEMPDIR = Path("DAJINResults", ".tempdir", NAME)
IS_CACHE_CONTROL = preprocess.validate_inputs.exists_cached_control(CONTROL, TEMPDIR)
IS_CACHE_GENOME = preprocess.validate_inputs.exists_cached_genome(GENOME, TEMPDIR, IS_CACHE_CONTROL)
UCSC_URL, GOLDENPATH_URL = None, None
if GENOME and not IS_CACHE_GENOME:
    UCSC_URL, GOLDENPATH_URL = preprocess.validate_inputs.check_and_fetch_genome(GENOME)

##########################################################
# Format inputs
##########################################################
SAMPLE_NAME = preprocess.format_inputs.extract_basename(SAMPLE)
CONTROL_NAME = preprocess.format_inputs.extract_basename(CONTROL)
FASTA_ALLELE = preprocess.format_inputs.dictionize_allele(ALLELE)

preprocess.format_inputs.make_directories(TEMPDIR)

if GENOME:
    GENOME_COODINATES, CHROME_SIZE = preprocess.format_inputs.get_coodinates_and_chromsize(
        TEMPDIR, GENOME, FASTA_ALLELE, UCSC_URL, GOLDENPATH_URL, IS_CACHE_GENOME
    )


flag1 = Path(TEMPDIR, "midsv", f"{SAMPLE_NAME}_control.jsonl").exists()
flag2 = Path(TEMPDIR, "midsv", f"{CONTROL_NAME}_control.jsonl").exists()
flag = flag1 and flag2

if not flag:
    ################################################################################
    # Export fasta files as single-FASTA format
    ################################################################################
    # TODO: use yeild, not export
    for identifier, sequence in FASTA_ALLELE.items():
        contents = "\n".join([">" + identifier, sequence]) + "\n"
        output_fasta = Path(TEMPDIR, "fasta", f"{identifier}.fasta")
        output_fasta.write_text(contents)
    ###############################################################################
    # Mapping with minimap2/mappy
    ###############################################################################
    for path_fasta in Path(TEMPDIR, "fasta").glob("*.fasta"):
        name_fasta = path_fasta.stem
        if name_fasta not in set(FASTA_ALLELE.keys()):
            continue
        if not IS_CACHE_CONTROL:
            preprocess.mappy_align.output_sam(TEMPDIR, path_fasta, name_fasta, CONTROL, CONTROL_NAME)
        preprocess.mappy_align.output_sam(TEMPDIR, path_fasta, name_fasta, SAMPLE, SAMPLE_NAME)
    ########################################################################
    # MIDSV conversion
    ########################################################################
    if not IS_CACHE_CONTROL:
        for path_sam in Path(TEMPDIR, "sam").glob(f"{CONTROL_NAME}*"):
            preprocess.call_midsv.output_midsv(TEMPDIR, path_sam, FASTA_ALLELE)
    for path_sam in Path(TEMPDIR, "sam").glob(f"{SAMPLE_NAME}*"):
        preprocess.call_midsv.output_midsv(TEMPDIR, path_sam, FASTA_ALLELE)
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
import midsv

# ?==============================================
dict_cssplit_control = defaultdict(list[dict])
for allele in FASTA_ALLELE.keys():
    path_control = Path(TEMPDIR, "midsv", f"{CONTROL_NAME}_{allele}.jsonl")
    cssplit_control = [cs["CSSPLIT"] for cs in midsv.read_jsonl(path_control)]
    dict_cssplit_control[allele] = cssplit_control

allele = "control"
sv = False

cssplit_sample = [cs["CSSPLIT"] for cs in classif_sample if cs["ALLELE"] == allele and cs["SV"] == sv]
cssplit_control = dict_cssplit_control[allele]
sequence = FASTA_ALLELE[allele]
masks_control = MASKS_CONTROL[allele]
diffloci, repetitive_delloci = screen_different_loci(
    cssplit_sample, cssplit_control, sequence, masks_control, alpha=0.01, threshold=0.01
)


def transpose(cssplit):
    cssplit = [cs.split(",") for cs in cssplit]
    return list(zip(*cssplit))


cssplit_sample_transposed = transpose(cssplit_sample)
cssplit_control_transposed = transpose(cssplit_control)

idx = 51
idx = 828
from collections import Counter
from collections import defaultdict

Counter(cssplit_sample_transposed[idx])
Counter(cssplit_control_transposed[idx])
x = len(cssplit_sample_transposed[idx])
y = len(cssplit_control_transposed[idx])
chistatistic([79, x - 1994 - 79], [200, y - 8296 - 200], threshold)
1994 / x * 100
8296 / y * 100


def count_percentage(cs):
    coverage = len(cs)
    d = defaultdict(float)
    for c in cs:
        d[c] += 1 / coverage * 100
    d = sorted(d.items(), key=lambda x: -x[1])
    return d


count_percentage(cssplit_sample_transposed[idx])
count_percentage(cssplit_control_transposed[idx])
Counter(cssplit_sample_transposed[idx])
Counter(cssplit_control_transposed[idx])

sample_cs, control_cs = cssplit_sample[idx]
alpha = 0.01
threshold = 0.01
table_sample, table_control = make_table(cssplit_sample, cssplit_control)
s, c = table_sample[idx], table_control[idx]
coverage = min(len(cssplit_control), len(cssplit_sample))
chistatistic([coverage, sum(s)], [coverage, sum(c)], threshold)
pval_del = chistatistic([coverage, s[0]], [coverage, c[0]], threshold)
pval_other = chistatistic([coverage, s[1]], [coverage, c[1]], threshold)
chistatistic([94, 2.3, 100 - 94 - 2.3], [91.4, 3.6, 100 - 91.4 - 3.6], threshold)

pval = pval_del < alpha or pval_other < alpha
# ?==============================================

MASKS_CONTROL = clustering.find_knockin_loci(TEMPDIR, FASTA_ALLELE, CONTROL_NAME)

DIFFLOCI_ALLELES, REPETITIVE_DELLOCI = clustering.extract_different_loci(
    TEMPDIR, classif_sample, MASKS_CONTROL, FASTA_ALLELE, CONTROL_NAME
)

clust_sample = clustering.add_labels(classif_sample, DIFFLOCI_ALLELES)

clust_sample = clustering.add_readnum(clust_sample)

clust_sample = clustering.add_percent(clust_sample)

clust_sample = clustering.update_labels(clust_sample)

########################################################################
# Consensus call
########################################################################

RESULT_SAMPLE, cons_percentage, cons_sequence = consensus.call(
    clust_sample, DIFFLOCI_ALLELES, REPETITIVE_DELLOCI, FASTA_ALLELE
)

d = defaultdict(int)
for res in RESULT_SAMPLE:
    d[res["NAME"]] += 1

d


# cons_percentage["allele3_albino_mutated_1.596%"][51]
# cssplit_sample = [cs for cs in clust_sample if cs["LABEL"] == 2]
# keys = []
# keys.append("control")
# keys.append(False)
# idx = 828

# ----------------------------------------------------------
# Conseusns Report：FASTA/HTML/VCF
# ----------------------------------------------------------
# FASTA
for header, cons_seq in cons_sequence.items():
    cons_fasta = postprocess.report_files.to_fasta(header, cons_seq)
    Path(TEMPDIR, "report", "FASTA", f"{SAMPLE_NAME}_{header}.fasta").write_text(cons_fasta)

# HTML
for header, cons_per in cons_percentage.items():
    cons_html = postprocess.report_files.to_html(header, cons_per)
    Path(TEMPDIR, "report", "HTML", f"{SAMPLE_NAME}_{header}.html").write_text(cons_html)

# BAM and igvjs
postprocess.report_bam.output_bam_control(TEMPDIR, CONTROL_NAME, GENOME, GENOME_COODINATES, CHROME_SIZE, THREADS)
postprocess.report_bam.output_bam_sample(
    TEMPDIR, RESULT_SAMPLE, SAMPLE_NAME, GENOME, GENOME_COODINATES, CHROME_SIZE, THREADS
)

# VCF
# working in progress


########################################################################
# MEMO
########################################################################


# ALLELE = "control"
# SV = False
# cssplit_sample = []
# for group in classif_sample:
#     if group["ALLELE"] == ALLELE and group["SV"] == SV:
#         cssplit_sample.append(group["CSSPLIT"])

# scores = [[1, 1, 1, 1], [1, 1, 1, 1]]
# set(scores)
# any(scores)
# np.all(scores)
