from itertools import groupby
from pathlib import Path

file_input = "examples/flox-cables2/AyabeTask1/batch.csv"

inputs = Path("examples/flox-cables2/AyabeTask1/batch.csv").read_text().strip().split("\n")

###############################################################################
# Check and format inputs (batch.csv)
###############################################################################

headers = inputs[0].split(",")

headers_argments = set(["sample", "control", "allele", "output", "genome", "debug", "threads"])
for header in headers:
    if header not in headers_argments:
        raise ValueError("Column names must be 'sample', 'control', 'allele', 'output', 'genome', and 'threads'")

index_control = [i for i, header in enumerate(headers) if header == "control"][0]
contents = [i.split(",") for i in inputs[1:]]
contents.sort(key=lambda x: x[index_control])
for _, groups in groupby(contents, key=lambda x: x[index_control]):
    for group in groups:
        arguments = {h: g for h, g in zip(headers, group)}
        sample = arguments["sample"]
        control = arguments["control"]
        allele = arguments["allele"]
        output = arguments["output"]
        try:
            genome = arguments["genome"]
        except KeyError:
            genome = ""
        try:
            debug = arguments["debug"]
        except KeyError:
            debug = False
        try:
            threads = arguments["threads"]
        except KeyError:
            threads = 1

