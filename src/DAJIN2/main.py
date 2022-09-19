from __future__ import annotations
from itertools import groupby
import shutil
from pathlib import Path

# Custom modules
from src.DAJIN2.utils import argparser


#! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# def main():

########################################################################
# Argument parse
########################################################################

sample, control, allele, output, genome, debug, threads = argparser.parse()

# #* Point mutation
# sample, control, allele, output, genome, debug, threads = (
#     "examples/pm-tyr/barcode31.fq.gz",
#     "examples/pm-tyr/barcode32.fq.gz",
#     "examples/pm-tyr/design_tyr.fa",
#     "DAJIN_results",
#     "mm10",
#     True,
#     14
#     )


# #* 2-cut deletion
# sample, control, allele, output, genome, debug, threads = (
#     "examples/del-stx2/barcode25.fq.gz",
#     "examples/del-stx2/barcode30.fq.gz",
#     "examples/del-stx2/design_stx2.fa",
#     "DAJIN_results",
#     "mm10",
#     True,
#     14,
# )

# #* flox insertion
# sample, control, allele, output, genome, debug, threads = (
#     "examples/flox-cables2/AyabeTask1/barcode31.fq.gz",
#     "examples/flox-cables2/AyabeTask1/barcode42.fq.gz",
#     "examples/flox-cables2/AyabeTask1/design_cables2.fa",
#     "DAJIN_results",
#     "mm10",
#     True,
#     14,
# )

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
# Check and format inputs (sample/control/allele)
###############################################################################

from pathlib import Path
from importlib import reload
from src.DAJIN2.preprocess import check_inputs
from src.DAJIN2.preprocess import format_inputs
from urllib.error import URLError

reload(check_inputs)
reload(format_inputs)

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
    if flag_fail:
        raise URLError("UCSC Download Servers are currently down")
    goldenpath_url, flag_fail = check_inputs.available_url(goldenpath_urls)
    # Check input genome
    check_inputs.available_genome(genome, ucsc_url)


# ------------------------------------------------------------------------------
# Formats
# ------------------------------------------------------------------------------

sample_name = format_inputs.extract_basename(sample)
control_name = format_inputs.extract_basename(control)
dict_allele = format_inputs.dictionize_allele(allele)

if genome:
    genome_coodinates = format_inputs.fetch_coodinate(genome, ucsc_url, dict_allele["control"])
    chrom_size = format_inputs.fetch_chrom_size(genome_coodinates["chr"], genome, goldenpath_url)

########################################################################
# Make directories
########################################################################

from pathlib import Path
import shutil

if Path(output).exists():
    shutil.rmtree(output)

Path(output).mkdir(exist_ok=True)

if Path(".tmpDAJIN").exists():
    shutil.rmtree(".tmpDAJIN")

subdirectoris = ["fasta", "fastq", "sam", "midsv", "bam", "reports"]
for subdir in subdirectoris:
    Path(".tmpDAJIN", subdir).mkdir(parents=True, exist_ok=True)

################################################################################
# Export fasta files as single-FASTA format
################################################################################

# TODO: use yeild, not export
from pathlib import Path

for header, sequence in dict_allele.items():
    contents = "\n".join([">" + header, sequence]) + "\n"
    output_fasta = Path(".tmpDAJIN", "fasta", f"{header}.fasta")
    output_fasta.write_text(contents)


###############################################################################
# Mapping with minimap2/mappy
###############################################################################

from pathlib import Path
from src.DAJIN2.preprocess import mappy_align
from importlib import reload

reload(mappy_align)

for input_fasta in Path(".tmpDAJIN", "fasta").glob("*.fasta"):
    fasta_name = input_fasta.name.replace(".fasta", "")
    for fastq, fastq_name in zip([control, sample], [control_name, sample_name]):
        # Todo: 並行処理で高速化！！
        SAM = mappy_align.to_sam(str(input_fasta), fastq)
        output_sam = Path(".tmpDAJIN", "sam", f"{fastq_name}_{fasta_name}.sam")
        output_sam.write_text("\n".join(SAM))

########################################################################
# MIDSV conversion
########################################################################

from pathlib import Path
import midsv

for sampath in Path(".tmpDAJIN", "sam").iterdir():
    output = Path(".tmpDAJIN", "midsv", f"{sampath.stem}.jsonl")
    sam = midsv.read_sam(sampath)
    midsv_jsonl = midsv.transform(sam, midsv=False, cssplit=True, qscore=False)
    midsv.write_jsonl(midsv_jsonl, output)


########################################################################
# Phread scoreが0.1以下のリードについて再分配する
########################################################################

# Under construction...

########################################################################
# Classify alleles
########################################################################

from src.DAJIN2 import classification
from importlib import reload

reload(classification)

classif_sample = classification.classify_alleles(sample_name)
classif_control = classification.classify_alleles(control_name)

########################################################################
# Detect Structural variants
########################################################################

from src.DAJIN2 import classification
from importlib import reload

reload(classification)

for classifs in [classif_sample, classif_control]:
    for classif in classifs:
        if classification.detect_sv(classif["CSSPLIT"], threshold=50):
            classif["SV"] = True
        else:
            classif["SV"] = False

########################################################################
# Clustering
########################################################################

# -----------------------------------------------------------------------
# Extract significantly different base loci between Sample and Control
# -----------------------------------------------------------------------

import midsv
from pathlib import Path
from collections import defaultdict
from importlib import reload
from itertools import groupby
from src.DAJIN2 import clustering

reload(clustering)

dict_cssplit_control = defaultdict(list[dict])
for ALLELE in dict_allele.keys():
    path_control = Path(".tmpDAJIN", "midsv", f"{control_name}_{ALLELE}.jsonl")
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

from src.DAJIN2 import clustering
from copy import deepcopy

labels = []
label_start = 1
for (ALLELE, SV), group in groupby(classif_sample, key=lambda x: (x["ALLELE"], x["SV"])):
    key = f'{{"ALLELE": "{ALLELE}", "SV": {SV}}}'
    cssplit_sample = [g["CSSPLIT"] for g in group]
    diffloci = diffloci_by_alleles[key]
    scores = list(clustering.make_scores(cssplit_sample, diffloci))
    labels += [label + label_start for label in clustering.clustering(scores).tolist()]
    label_start = len(set(labels)) + 1

clust_sample = deepcopy(classif_sample)
for clust, label in zip(clust_sample, labels):
    clust["LABEL"] = label
    del clust["CSSPLIT"]

clust_sample.sort(key=lambda x: (x["ALLELE"], x["SV"], x["LABEL"]))

########################################################################
# Consensus call
########################################################################

from src.DAJIN2.consensus import module_consensus as consensus
from collections import defaultdict
from importlib import reload

reload(consensus)

path = Path(".tmpDAJIN", "midsv", f"{control_name}_control.jsonl")
cssplit_control = midsv.read_jsonl(path)

path = Path(".tmpDAJIN", "midsv", f"{sample_name}_control.jsonl")
cssplit_sample = midsv.read_jsonl(path)
cssplit_sample = consensus.join_listdicts(clust_sample, cssplit_sample, key="QNAME")
cssplit_sample.sort(key=lambda x: (x["ALLELE"], x["SV"], x["LABEL"]))

cons_percentage = defaultdict(list)
cons_sequence = defaultdict(list)
for (ALLELE, SV, LABEL), cssplits in groupby(cssplit_sample, key=lambda x: (x["ALLELE"], x["SV"], x["LABEL"])):
    cons_per = consensus.call_percentage(list(cssplits), cssplit_control)
    cons_seq = consensus.call_sequence(cons_per)
    key = f'{{"ALLELE": "{ALLELE}", "SV": {SV}, "LABEL": {LABEL}}}'
    cons_percentage[key] = cons_per
    cons_sequence[key] = cons_seq

########################################################################
# Report：アレル割合
# sample, allele name, #read, %read
########################################################################
import warnings

warnings.simplefilter("ignore")
from collections import defaultdict
from src.DAJIN2.report import report_af

reload(report_af)

clust_sample = report_af.call_allele_name(clust_sample, cons_sequence, dict_allele)

# ----------------------------------------------------------
# All data
# ----------------------------------------------------------

df_clust_sample = report_af.all_allele(clust_sample, sample_name)
df_clust_sample.to_csv(".tmpDAJIN/reports/read_classification_all.csv", index=False)
df_clust_sample.to_excel(".tmpDAJIN/reports/read_classification_all.xlsx", index=False)

# ----------------------------------------------------------
# Summary data
# ----------------------------------------------------------

df_allele_frequency = report_af.summary_allele(clust_sample, sample_name)
df_allele_frequency.to_csv(".tmpDAJIN/reports/read_classification_summary.csv", index=False)
df_allele_frequency.to_excel(".tmpDAJIN/reports/read_classification_summary.xlsx", index=False)

# ----------------------------------------------------------
# Visualization
# ----------------------------------------------------------
g = report_af.plot(df_allele_frequency)
g.save(filename=".tmpDAJIN/reports/tmp_output.png", dpi=350)
g.save(filename=".tmpDAJIN/reports/tmp_output.pdf")

########################################################################
# Report：各種ファイルフォーマットに出力
########################################################################

from src.DAJIN2.report import report_files

reload(report_files)

# FASTA
headers = df_allele_frequency["allele name"].to_list()
for heaer, cons_seq in zip(headers, cons_sequence.values()):
    cons_fasta = report_files.to_fasta(heaer, cons_seq)
    Path(f".tmpDAJIN/reports/{heaer}.fasta").write_text(cons_fasta)

# HTML
for heaer, cons_per in zip(headers, cons_percentage.values()):
    cons_html = report_files.to_html(heaer, cons_per)
    Path(f".tmpDAJIN/reports/{heaer}.html").write_text(cons_html)


# VCF


########################################################################
# Report：BAM
########################################################################

import pysam
from collections import defaultdict
from src.DAJIN2.report import report_bam

reload(report_bam)

Path(".tmpDAJIN", "bam", "tmp_sam").mkdir(parents=True, exist_ok=True)

path_sam = Path(".tmpDAJIN", "sam", f"{sample_name}_control.sam")
sam = midsv.read_sam(path_sam)

sam_update = sam.copy()
sam_update = report_bam.remove_overlapped_reads(sam_update)
sam_update = report_bam.remove_microhomology(sam_update)

path_sam = f".tmpDAJIN/bam/tmp_sam/{sample_name}_control_updated.sam"
if genome:
    sam_update = report_bam.realign(sam_update, genome_coodinates, chrom_size)
    report_bam.write_sam(sam_update, path_sam)
else:
    report_bam.write_sam(sam_update, path_sam)

pysam.sort("-@", f"{threads}", "-o", f".tmpDAJIN/bam/{sample_name}.bam", path_sam)
pysam.index("-@", f"{threads}", f".tmpDAJIN/bam/{sample_name}.bam")

sam_headers = [s for s in sam_update if s[0].startswith("@")]
sam_contents = [s for s in sam_update if not s[0].startswith("@")]
sam_groups = report_bam.group_by_name(sam_contents, clust_sample)
for key, sam_content in sam_groups.items():
    path_sam = f".tmpDAJIN/bam/tmp_sam/{key}.sam"
    report_bam.write_sam(sam_headers + sam_content, path_sam)
    pysam.sort("-@", f"{threads}", "-o", f".tmpDAJIN/bam/{sample_name}_{key}.bam", path_sam)
    pysam.index("-@", f"{threads}", f".tmpDAJIN/bam/{sample_name}_{key}.bam")


########################################################################
# Report：IGV.js
########################################################################

# IGV.js（10本のリードのみ表示）のために各サンプル50本程度のリードのみを抽出する

for path_sam in Path(".tmpDAJIN", "sam").glob("*_control.sam"):
    sam = midsv.read_sam(path_sam)
    pysam.sort("-o", f".tmpDAJIN/bam/{name}.bam", str(path_sam))
    pysam.index(f".tmpDAJIN/bam/{name}.bam")


########################################################################
# Finish call
########################################################################

if not debug:
    shutil.rmtree(".tmpDAJIN")

# if __name__ == "__main__":
#     main()
