from __future__ import annotations

import shutil
from pathlib import Path

import midsv
import pandas as pd
import plotly.express as px

from DAJIN2.utils.config import DAJIN_RESULTS_DIR, TEMP_ROOT_DIR


def extract_all_info(batch_directory: Path) -> pd.DataFrame:
    df_results = pd.DataFrame()
    colum = ["SAMPLE", "QNAME", "NAME", "READNUM", "PERCENT", "LABEL", "ALLELE"]
    for path_result in sorted(batch_directory.iterdir()):
        sample_name = path_result.stem
        result_jsonl = midsv.read_jsonl(path_result)
        df_clust_sample = pd.DataFrame(result_jsonl)
        df_clust_sample["SAMPLE"] = sample_name
        df_clust_sample = df_clust_sample[colum]
        df_results = pd.concat([df_results, df_clust_sample])
    return df_results


def summarize_info(df_all: pd.DataFrame) -> pd.DataFrame:
    df_summary = df_all.drop(columns=["QNAME", "LABEL", "ALLELE"]).drop_duplicates()
    return df_summary


def output_plot(df_summary: pd.DataFrame, report_directory: Path):
    df_plot = df_summary.copy()
    name = df_plot.NAME.tolist()
    names = [n.split("_") for n in name]
    alleletype = [n[1] + " " + n[2] for n in names]
    df_plot["ALLELETYPE"] = alleletype
    fig = px.bar(
        df_plot,
        x="SAMPLE",
        y="PERCENT",
        color="ALLELETYPE",
        text="PERCENT",
        labels={"SAMPLE": "Samples", "PERCENT": "% of reads", "ALLELETYPE": "Alelle type"},
    )
    fig.update_traces(textposition="inside", cliponaxis=False)
    fig.update_xaxes(categoryorder="category ascending")
    output_filename = str(Path(report_directory, "read_plot"))
    fig.write_html(f"{output_filename}.html")
    # if kaleido is installed, output a pdf
    try:
        fig.write_image(f"{output_filename}.pdf")
    except ValueError:
        pass


def report(name: str):
    report_directory = Path(DAJIN_RESULTS_DIR, name)
    report_directory.mkdir(exist_ok=True, parents=True)
    shutil.copytree(Path(TEMP_ROOT_DIR, name, "report"), report_directory, dirs_exist_ok=True)
    df_all = extract_all_info(Path(TEMP_ROOT_DIR, name, "result"))
    df_all.to_csv(Path(report_directory, "read_all.csv"), index=False)
    df_summary = summarize_info(df_all)
    df_summary.to_csv(Path(report_directory, "read_summary.csv"), index=False)
    output_plot(df_summary, report_directory)
