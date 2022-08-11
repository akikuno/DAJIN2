# * flox insertion
sample, control, allele, output, genome, debug, threads = (
    "examples/flox-cables2/AyabeTask1/barcode31.fq.gz",
    "examples/flox-cables2/AyabeTask1/barcode42.fq.gz",
    "examples/flox-cables2/AyabeTask1/design_cables2.fa",
    "DAJIN_results",
    "mm10",
    True,
    14,
)

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


classif_sample.sort(key=lambda x: (x["ALLELE"], x["SV"]))
classif_sample_groupby = groupby(classif_sample, key=lambda x: (x["ALLELE"], x["SV"]))

allele_diffloci = defaultdict(list[dict])
dict_control_cssplit = defaultdict(list[dict])
for (ALLELE, SV), group in classif_sample_groupby:
    sample_cssplit = [record["CSSPLIT"] for record in group]
    if not dict_control_cssplit[ALLELE]:
        control_path = Path(".tmpDAJIN", "midsv", f"{control_name}_{ALLELE}.jsonl")
        control_cssplit = [cs["CSSPLIT"] for cs in midsv.read_jsonl(control_path)]
        dict_control_cssplit[ALLELE] = control_cssplit
    control_cssplit = dict_control_cssplit[ALLELE]
    sequence = dict_allele[ALLELE]
    diffloci = clustering.screen_different_loci(sample_cssplit, control_cssplit, sequence, alpha=0.01, threshold=0.05)
    allele_diffloci[f'{{"ALLELE": "{ALLELE}", "SV": {SV}}}'] = diffloci


# -----------------------------------------------------------------------
# Annotate scores to sample's reads
# -----------------------------------------------------------------------

from src.DAJIN2 import clustering
from collections import defaultdict

reload(clustering)

sample_scores_by_allele = defaultdict(list[dict])
for (ALLELE, SV), group in groupby(classif_sample, key=lambda x: (x["ALLELE"], x["SV"])):
    key = f'{{"ALLELE": "{ALLELE}", "SV": {SV}}}'
    sample_cssplit = [g["CSSPLIT"] for g in group]
    diff_loci = allele_diffloci[key]
    sample_scores_by_allele[key] = clustering.make_scores(sample_cssplit, diff_loci)

