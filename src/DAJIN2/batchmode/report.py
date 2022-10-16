from __future__ import annotations
import pandas as pd
import plotly.express as px
from plotnine import ggplot, aes, geom_bar, theme, theme_bw, element_blank, labs, scale_y_continuous
import midsv
from pathlib import Path


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


# def plot(df_summary: pd.DataFrame):
#     df_plot = df_summary.copy()
#     name = df_plot.NAME.tolist()
#     names = [n.split("_") for n in name]
#     alleletype = [n[1] + " " + n[2] for n in names]
#     df_plot["ALLELETYPE"] = alleletype
#     g = (
#         ggplot(df_plot, aes(x="SAMPLE", y="READNUM", fill="ALLELETYPE"))
#         + geom_bar(position="fill", stat="identity", colour="black")
#         + scale_y_continuous(labels=lambda l: ["%d%%" % (v * 100) for v in l])
#         + theme_bw()
#         + theme(axis_title_x=element_blank())
#         + labs(fill="Allele type", y="Percentage of reads")
#     )
#     return g


###############################################################################
# igvjs
###############################################################################

