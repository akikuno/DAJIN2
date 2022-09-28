from __future__ import annotations
import warnings
from itertools import groupby
from pathlib import Path
from src.DAJIN2.core import main
from importlib import reload

warnings.simplefilter("ignore")
reload(main)

###############################################################################
# Parse arguments (batch.csv)
###############################################################################

input_file = "examples/flox-cables2/AyabeTask1/batch.csv"

inputs = Path(input_file).read_text().strip().split("\n")

keys = inputs[0].split(",")

keys_set = set(keys)
keys_required = set(["sample", "control", "allele"])
keys_all = set(["sample", "control", "allele", "output", "genome", "threads", "debug"])

if not keys_required.issubset(keys_set):
    raise ValueError(f"Column names of {input_file} must contain 'sample', 'control', and 'allele'")

if not keys_set.issubset(keys_all):
    raise ValueError(
        f"Accepted column names of {input_file} are 'sample', 'control', 'allele', 'output', 'genome', and 'threads'."
    )

index_control = [i for i, header in enumerate(keys) if header == "control"][0]
contents = [i.split(",") for i in inputs[1:]]
contents.sort(key=lambda x: x[index_control])
for _, groups in groupby(contents, key=lambda x: x[index_control]):
    for group in groups:
        arguments = {h: g for h, g in zip(keys, group)}
        # results_main = main(arguments)
        # main(arguments)


results_main = main.main(arguments)

