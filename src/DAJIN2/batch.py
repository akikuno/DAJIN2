from __future__ import annotations
from pathlib import Path
from itertools import groupby
import shutil
import pandas as pd
from .core import main
import wslPath
import midsv
from src.DAJIN2.batchmode import report
from jinja2 import Template, Environment, FileSystemLoader


def batch(arguments: dict[str]):
    path_batchfile = arguments["file"]
    threads = arguments["threads"]
    debug = arguments["debug"]

    ###############################################################################
    # Check arguments
    ###############################################################################
    try:
        path_batchfile = wslPath.toPosix(path_batchfile)
    except ValueError:
        pass

    # ----------------------------------------------------------
    # Check file exists
    # ----------------------------------------------------------
    if not Path(path_batchfile).exists():
        raise FileNotFoundError(f"'{path_batchfile}' does not exist.")

    try:
        df_batchfile = pd.read_excel(path_batchfile)
        inputs = []
        inputs = df_batchfile.columns.to_list()
        inputs += df_batchfile.values.tolist()
    except ValueError:
        inputs = [s.split(",") for s in Path(path_batchfile).read_text().strip().split("\n")]

    # ----------------------------------------------------------
    # Check Column
    # ----------------------------------------------------------
    keys = inputs[0]
    keys_set = set(keys)
    keys_required = set(["sample", "control", "allele", "name"])
    keys_all = set(["sample", "control", "allele", "name", "genome"])

    if not keys_required.issubset(keys_set):
        raise ValueError(f"{path_batchfile} must contain 'sample', 'control', and 'allele' in the header")

    if not keys_set.issubset(keys_all):
        raise ValueError(
            f"Accepted header names of {path_batchfile} are 'sample', 'control', 'allele', 'name', and 'genome'."
        )

    ##############################################################################
    # Perform DAJIN
    ##############################################################################
    index_name = [i for i, key in enumerate(keys) if key == "name"][0]
    contents = inputs[1:]
    contents.sort(key=lambda x: x[index_name])

    for name, groups in groupby(contents, key=lambda x: x[index_name]):
        report_directory = Path("DAJINResults", name)
        report_directory.mkdir(exist_ok=True, parents=True)
        batch_directory = Path("DAJINResults", ".tempdir", name, "batch")
        batch_directory.mkdir(exist_ok=True, parents=True)
        for group in groups:
            args = {h: g for h, g in zip(keys, group)}
            args["threads"] = threads
            # perform DAJIN
            SAMPLE_NAME, RESULT = main.main(args)
            # output the result
            midsv.write_jsonl(RESULT, Path(batch_directory, f"{SAMPLE_NAME}.jsonl"))
            # ? DEVELOPMENT------------------------------------------------------------------
            print(f"DEVELOPLENT: batch.py: {SAMPLE_NAME}.jsonl is completed")
            # ? DEVELOPMENT------------------------------------------------------------------
        # Reports -----------------------------------------
        for dir in Path("DAJINResults", ".tempdir", name, "report").iterdir():
            shutil.copytree(dir, Path(report_directory, dir.stem), dirs_exist_ok=True)
        df_all = report.all_info(batch_directory)
        df_all.to_csv(Path(report_directory, "read_all.csv"), index=False)
        df_summary = report.summary_info(df_all)
        df_summary.to_csv(Path(report_directory, "read_summary.csv"), index=False)
        report.output_plot(df_summary, report_directory)


##############################################################################
# Reports
##############################################################################

# IGV.js

# df_summary = df_results.drop(columns=["QNAME", "LABEL", "ALLELE", "SV"]).drop_duplicates()
# env = Environment(loader=FileSystemLoader("./", encoding="utf8"))
# template = env.get_template("src/DAJIN2/template_igvjs.html")  #!? DEVELOPMENT ==================================
# # template = env.get_template("template_igvjs.html")
# params_genome = dict()
# params_reference = dict()
# params_tracks = dict()
# params_genome = {"genome": {"exist": False}}
# if GENOME:
#     chrom, start, end, _ = GENOME_COODINATES.values()
#     params_genome = {"genome": {"exist": True, "genome": GENOME, "locus": f"{chrom}:{start}-{end}"}}
# else:
#     params_reference = {"reference": {"urlfa": f""}}

# params_reference = {"reference": {"urlfa": "@@@@@@@@@@@@@@@@@@@"}}

# params_tracks = """{
#     "tracks":
#         {{ contents }},
# }"""
# params_tracks = Template(params_tracks)

# SAMPLENAMES = ["hoge", "fuga"]
# URLBAM = ["DAJINResults/.tempdir/Ayabe-Task1/bam/barcode36_allele1_flox_mutated_50.3%.bam", "fuga2"]
# URLBAI = ["foooo", "foooo2"]
# contents = [{"samplename": s, "urlbam": b, "urlbai": i} for s, b, i in zip(SAMPLENAMES, URLBAM, URLBAI)]

# params_tracks = eval(params_tracks.render(contents=contents))

# params_genome.update(**params_reference, **params_tracks)
# params = params_genome.copy()

# result = template.render(params)
# print(result)