from __future__ import annotations
import pandas as pd
import shutil
from pathlib import Path
import midsv
import plotly.express as px


def all_info(batch_directory: Path) -> pd.DataFrame:
    df_results = pd.DataFrame()
    colum = ["SAMPLE", "QNAME", "NAME", "READNUM", "PERCENT", "LABEL", "ALLELE", "SV"]
    for path_result in batch_directory.iterdir():
        sample_name = path_result.stem
        result_jsonl = midsv.read_jsonl(path_result)
        df_clust_sample = pd.DataFrame(result_jsonl)
        df_clust_sample["SAMPLE"] = sample_name
        df_clust_sample = df_clust_sample[colum]
        df_results = pd.concat([df_results, df_clust_sample])
    return df_results


def summary_info(df_all: pd.DataFrame) -> pd.DataFrame:
    df_summary = df_all.drop(columns=["QNAME", "LABEL", "ALLELE", "SV"]).drop_duplicates()
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
    output_filename = str(Path(report_directory, "read_plot"))
    fig.write_html(f"{output_filename}.html")
    # if kaleido is installed, output a pdf
    try:
        fig.write_image(f"{output_filename}.pdf")
    except ValueError:
        pass


def report(name: str):
    report_directory = Path("DAJINResults", name)
    report_directory.mkdir(exist_ok=True, parents=True)
    for dir in Path("DAJINResults", ".tempdir", name, "report").iterdir():
        shutil.copytree(dir, Path(report_directory, dir.stem), dirs_exist_ok=True)
    df_all = all_info(Path("DAJINResults", ".tempdir", name, "result"))
    df_all.to_csv(Path(report_directory, "read_all.csv"), index=False)
    df_summary = summary_info(df_all)
    df_summary.to_csv(Path(report_directory, "read_summary.csv"), index=False)
    output_plot(df_summary, report_directory)
