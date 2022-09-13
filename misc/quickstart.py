import midsv
from pathlib import Path
import shutil
from src.DAJIN2.preprocess import mappy_align
from importlib import reload

reload(mappy_align)

# * Point mutation
sample, control, allele, output, genome, debug, threads = (
    "examples/pm-tyr/barcode31.fq.gz",
    "examples/pm-tyr/barcode32.fq.gz",
    "examples/pm-tyr/design_tyr.fa",
    "DAJIN_results",
    "mm10",
    True,
    14,
)


# # * 2-cut deletion
# sample, control, allele, output, genome, debug, threads = (
#     "examples/del-stx2/barcode25.fq.gz",
#     "examples/del-stx2/barcode30.fq.gz",
#     "examples/del-stx2/design_stx2.fa",
#     "DAJIN_results",
#     "mm10",
#     True,
#     14,
# )

# # * flox insertion
# sample, control, allele, output, genome, debug, threads = (
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
from src.DAJIN2.preprocess import format_input

reload(format_input)

# ------------------------------------------------------------------------------
# Check input path
# ------------------------------------------------------------------------------

if not Path(control).exists():
    raise FileNotFoundError(f"{control} is not found")
elif not Path(sample).exists():
    raise FileNotFoundError(f"{sample} is not found")
elif not Path(allele).exists():
    raise FileNotFoundError(f"{allele} is not found")

# ------------------------------------------------------------------------------
# Check formats (extensions and contents)
# ------------------------------------------------------------------------------

for file_name in sample, control:
    format_input.check_fastq_extension(file_name)

for file_name in sample, control:
    format_input.check_fastq_content(file_name)

format_input.check_fasta_content(allele)

# ------------------------------------------------------------------------------
# Formats
# ------------------------------------------------------------------------------

sample_name = format_input.extract_basename(sample)
control_name = format_input.extract_basename(control)
dict_allele = format_input.dictionize_allele(allele)

flag1 = Path(".tmpDAJIN", "midsv", f"{sample_name}_control.jsonl").exists()
flag2 = Path(".tmpDAJIN", "midsv", f"{control_name}_control.jsonl").exists()
flag = flag1 and flag2

if not flag:
    if Path(output).exists():
        shutil.rmtree(output)
    Path(output).mkdir(exist_ok=True)
    if Path(".tmpDAJIN").exists():
        shutil.rmtree(".tmpDAJIN")
    for subdir in ["fasta", "fastq", "sam", "midsv"]:
        Path(".tmpDAJIN", subdir).mkdir(parents=True, exist_ok=True)
    # TODO: use yeild, not export
    for header, sequence in dict_allele.items():
        contents = "\n".join([">" + header, sequence]) + "\n"
        output_fasta = Path(".tmpDAJIN", "fasta", f"{header}.fasta")
        output_fasta.write_text(contents)
    ###############################################################################
    # Mapping with minimap2/mappy
    ###############################################################################
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
    for sampath in Path(".tmpDAJIN", "sam").iterdir():
        output = Path(".tmpDAJIN", "midsv", f"{sampath.stem}.jsonl")
        sam = midsv.read_sam(sampath)
        midsv_jsonl = midsv.transform(sam, midsv=False, cssplit=True, qscore=False)
        midsv.write_jsonl(midsv_jsonl, output)


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

classif_sample.sort(key=lambda x: (x["ALLELE"], x["SV"]))

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


ALLELE, SV = "control", False
cssplit_sample = []
for c in classif_sample:
    if c["ALLELE"] == ALLELE and c["SV"] == SV:
        cssplit_sample.append(c["CSSPLIT"])

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

# -----------------------------------------------------------------------
# Memo: Clustering
# -----------------------------------------------------------------------

# from importlib import reload
# from collections import defaultdict
# from pprint import pprint

# d = defaultdict(int)
# for c in classif_sample:
#     ALLELE, SV = c["ALLELE"], c["SV"]
#     key = f'{{"ALLELE": "{ALLELE}", "SV": {SV}}}'
#     d[key] += 1

# pprint(d)

# # ? read check
# ALLELE = "control"
# SV = False
# cssplit_sample = []
# cssplit_control = dict_cssplit_control[ALLELE]
# qnames = []
# diffloci = diffloci_by_alleles[f'{{"ALLELE": "{ALLELE}", "SV": {SV}}}']
# for c in classif_sample:
#     if c["ALLELE"] == ALLELE and c["SV"] == SV:
#         cssplit_sample.append(c["CSSPLIT"])
#         qnames.append(c["QNAME"])


# scores = list(clustering.make_scores(cssplit_sample, cssplit_control, diffloci))

# labels = clustering.clustering(scores)

# # ? check

# from collections import Counter

# Counter(labels)

# scores[10]
# labels[11]
# scores[12]

# qa = []
# for i, (q, a) in enumerate(zip(qnames, labels.tolist())):
#     qa.append({"QNAME": q, "LABEL": a})

# from itertools import groupby

# qa.sort(key=lambda x: x["LABEL"])
# for LABEL, group in groupby(qa, key=lambda x: x["LABEL"]):
#     with open(f"tmp_{LABEL}", "a") as f:
#         f.write("^@\n")
#         for g in group:
#             # print(g["QNAME"])
#             qname = g["QNAME"]
#             f.write(f"{qname}\n")
