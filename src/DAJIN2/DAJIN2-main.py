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

sample, control, allele, output, genome, debug, threads = argparser.parse()

########################################################################
# Whether existing cached control
########################################################################

# CACHEDIR = os.path.join(tempfile.gettempdir(), "DAJIN")
# os.makedirs(CACHEDIR, exist_ok=True)

# if cache_control.exists(control, CACHEDIR):
#     IS_CACHED = True
# else:
#     cache_control.save_header(control, CACHEDIR)
#     IS_CACHED = False

###############################################################################
# Check inputs (sample/control/allele/genome)
###############################################################################

# ------------------------------------------------------------------------------
# Check input path
# ------------------------------------------------------------------------------

check_inputs.exists(control)
check_inputs.exists(sample)
check_inputs.exists(allele)

# ------------------------------------------------------------------------------
# Check formats (extensions and contents)
# ------------------------------------------------------------------------------

check_inputs.fastq_extension(control)
check_inputs.fastq_content(control)
check_inputs.fastq_extension(sample)
check_inputs.fastq_content(sample)
check_inputs.fasta_content(allele)

# ------------------------------------------------------------------------------
# Check genomes if genome is inputted
# ------------------------------------------------------------------------------

if genome:
    # Check UCSC Server
    ucsc_urls = [
        "https://genome.ucsc.edu/",
        "https://genome-asia.ucsc.edu/",
        "https://genome-euro.ucsc.edu/",
    ]
    ucsc_url, flag_fail = check_inputs.available_url(ucsc_urls)
    if flag_fail:
        raise URLError("UCSC Servers are currently down")
    # Check UCSC Download Server
    goldenpath_urls = [
        "https://hgdownload.cse.ucsc.edu/goldenPath",
        "http://hgdownload-euro.soe.ucsc.edu/goldenPath",
    ]
    goldenpath_url, flag_fail = check_inputs.available_url(goldenpath_urls)
    if flag_fail:
        raise URLError("UCSC Download Servers are currently down")
    # Check input genome
    check_inputs.available_genome(genome, ucsc_url)


########################################################################
# Make directories
########################################################################

if Path(output).exists():
    shutil.rmtree(output)

Path(output).mkdir(exist_ok=True)

if Path(output, ".tempdir").exists():
    shutil.rmtree(Path(output, ".tempdir"))

subdirectoris = ["fasta", "sam", "midsv", "bam", "reports"]
for subdir in subdirectoris:
    Path(output, ".tempdir", subdir).mkdir(parents=True, exist_ok=True)

###############################################################################
# Format inputs (sample/control/allele/genome)
###############################################################################

sample_name = format_inputs.extract_basename(sample)
control_name = format_inputs.extract_basename(control)
dict_allele = format_inputs.dictionize_allele(allele)

if genome:
    genome_coodinates = format_inputs.fetch_coodinate(genome, ucsc_url, dict_allele["control"])
    chrom_size = format_inputs.fetch_chrom_size(genome_coodinates["chr"], genome, goldenpath_url)

################################################################################
# Export fasta files as single-FASTA format
################################################################################

# TODO: use yeild, not export

for header, sequence in dict_allele.items():
    contents = "\n".join([">" + header, sequence]) + "\n"
    output_fasta = Path(output, ".tempdir", "fasta", f"{header}.fasta")
    output_fasta.write_text(contents)


###############################################################################
# Mapping with minimap2/mappy
###############################################################################


for input_fasta in Path(output, ".tempdir", "fasta").glob("*.fasta"):
    fasta_name = input_fasta.name.replace(".fasta", "")
    for fastq, fastq_name in zip([control, sample], [control_name, sample_name]):
        # Todo: 並行処理で高速化！！
        SAM = mappy_align.to_sam(str(input_fasta), fastq)
        output_sam = Path(output, ".tempdir", "sam", f"{fastq_name}_{fasta_name}.sam")
        output_sam.write_text("\n".join(SAM))

########################################################################
# MIDSV conversion
########################################################################

for sampath in Path(output, ".tempdir", "sam").iterdir():
    output_jsonl = Path(output, ".tempdir", "midsv", f"{sampath.stem}.jsonl")
    sam = midsv.read_sam(sampath)
    midsv_jsonl = midsv.transform(sam, midsv=False, cssplit=True, qscore=False)
    midsv.write_jsonl(midsv_jsonl, output_jsonl)

########################################################################
# Classify alleles
########################################################################


path_midsv = Path(output, ".tempdir", "midsv").glob(f"{sample_name}*")
classif_sample = classification.classify_alleles(path_midsv, sample_name)

path_midsv = Path(output, ".tempdir", "midsv").glob(f"{control_name}*")
classif_control = classification.classify_alleles(path_midsv, control_name)

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
for ALLELE in dict_allele.keys():
    path_control = Path(output, ".tempdir", "midsv", f"{control_name}_{ALLELE}.jsonl")
    cssplit_control = [cs["CSSPLIT"] for cs in midsv.read_jsonl(path_control)]
    dict_cssplit_control[ALLELE] = cssplit_control

classif_sample.sort(key=lambda x: (x["ALLELE"], x["SV"]))
diffloci_by_alleles = defaultdict(list[dict])
for (ALLELE, SV), group in groupby(classif_sample, key=lambda x: (x["ALLELE"], x["SV"])):
    cssplit_sample = [record["CSSPLIT"] for record in group]
    cssplit_control = dict_cssplit_control[ALLELE]
    sequence = dict_allele[ALLELE]
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


path = Path(output, ".tempdir", "midsv", f"{control_name}_control.jsonl")
cssplit_control = midsv.read_jsonl(path)

path = Path(output, ".tempdir", "midsv", f"{sample_name}_control.jsonl")
cssplit_sample = midsv.read_jsonl(path)
cssplit_sample = consensus.join_listdicts(clust_sample, cssplit_sample, key="QNAME")
cssplit_sample.sort(key=lambda x: (x["ALLELE"], x["SV"], x["LABEL"]))

cons_percentage = defaultdict(list)
cons_sequence = defaultdict(list)
for keys, cssplits in groupby(cssplit_sample, key=lambda x: (x["ALLELE"], x["SV"], x["LABEL"])):
    cssplits = list(cssplits)
    cons_per = consensus.call_percentage(cssplits, cssplit_control)
    cons_seq = consensus.call_sequence(cons_per)
    allele_name = consensus.call_allele_name(keys, cons_seq, dict_allele)
    cons_percentage[allele_name] = cons_per
    cons_sequence[allele_name] = cons_seq
    for cs in cssplits:
        cs["NAME"] = allele_name

for cs in cssplit_sample:
    del cs["RNAME"]
    del cs["CSSPLIT"]


# ----------------------------------------------------------
# Conseusns Report：FASTA/HTML/VCF
# ----------------------------------------------------------

# FASTA
for header, cons_seq in cons_sequence.items():
    cons_fasta = consensus.to_fasta(header, cons_seq)
    Path(f"{output}/.tempdir/reports/{sample_name}_{header}.fasta").write_text(cons_fasta)

# HTML
for header, cons_per in cons_percentage.items():
    cons_html = consensus.to_html(header, cons_per)
    Path(f"{output}/.tempdir/reports/{sample_name}_{header}.html").write_text(cons_html)


# VCF
# working in progress

########################################################################
# Report：アレル割合
# sample, allele name, #read, %read
########################################################################

# ----------------------------------------------------------
# All data
# ----------------------------------------------------------

df_clust_sample = report_af.all_allele(clust_sample, sample_name)
df_clust_sample.to_csv(f"{output}/.tempdir/reports/{sample_name}_read_classification_all.csv", index=False)
df_clust_sample.to_excel(f"{output}/.tempdir/reports/{sample_name}_read_classification_all.xlsx", index=False)

# ----------------------------------------------------------
# Summary data
# ----------------------------------------------------------

df_allele_frequency = report_af.summary_allele(clust_sample, sample_name)
df_allele_frequency.to_csv(f"{output}/.tempdir/reports/{sample_name}_read_classification_summary.csv", index=False)
df_allele_frequency.to_excel(f"{output}/.tempdir/reports/{sample_name}_read_classification_summary.xlsx", index=False)

# ----------------------------------------------------------
# Visualization
# ----------------------------------------------------------
g = report_af.plot(df_allele_frequency)
g.save(filename=f"{output}/.tempdir/reports/tmp_output.png", dpi=350)
g.save(filename=f"{output}/.tempdir/reports/tmp_output.pdf")


########################################################################
# Report：BAM
########################################################################


Path(output, ".tempdir", "bam", "tmp_sam").mkdir(parents=True, exist_ok=True)

path_sam = Path(output, ".tempdir", "sam", f"{sample_name}_control.sam")
sam = midsv.read_sam(path_sam)

sam_update = sam.copy()
sam_update = report_bam.remove_overlapped_reads(sam_update)
sam_update = report_bam.remove_microhomology(sam_update)

path_sam = f"{output}/.tempdir/bam/tmp_sam/{sample_name}_control_updated.sam"
if genome:
    sam_update = report_bam.realign(sam_update, genome_coodinates, chrom_size)
    report_bam.write_sam(sam_update, path_sam)
else:
    report_bam.write_sam(sam_update, path_sam)

pysam.sort("-@", f"{threads}", "-o", f"{output}/.tempdir/bam/{sample_name}.bam", path_sam)
pysam.index("-@", f"{threads}", f"{output}/.tempdir/bam/{sample_name}.bam")

sam_headers = [s for s in sam_update if s[0].startswith("@")]
sam_contents = [s for s in sam_update if not s[0].startswith("@")]
sam_groups = report_bam.group_by_name(sam_contents, clust_sample)
for key, sam_content in sam_groups.items():
    path_sam = f"{output}/.tempdir/bam/tmp_sam/{key}.sam"
    report_bam.write_sam(sam_headers + sam_content, path_sam)
    pysam.sort("-@", f"{threads}", "-o", f"{output}/.tempdir/bam/{sample_name}_{key}.bam", path_sam)
    pysam.index("-@", f"{threads}", f"{output}/.tempdir/bam/{sample_name}_{key}.bam")


########################################################################
# Report：IGV.js
########################################################################

# IGV.js（10本のリードのみ表示）のために各サンプル50本程度のリードのみを抽出する

# for path_sam in Path(output, ".tempdir", "sam").glob("*_control.sam"):
#     sam = midsv.read_sam(path_sam)
#     pysam.sort("-o", f"{output}/.tempdir/bam/{name}.bam", str(path_sam))
#     pysam.index(f"{output}/.tempdir/bam/{name}.bam")


########################################################################
# Finish call
########################################################################

if not debug:
    shutil.rmtree(Path(output, ".tempdir"))

# if __name__ == "__main__":
#     main()
