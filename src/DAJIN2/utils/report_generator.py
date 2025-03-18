from __future__ import annotations

import shutil
from pathlib import Path

import plotly.express as px

from DAJIN2.utils import io
from DAJIN2.utils.config import DAJIN_RESULTS_DIR, TEMP_ROOT_DIR


def rename_allele(allele: str) -> str:
    if allele == "control" or allele.startswith("insertion"):
        return allele.capitalize()
    return allele


def rename_type(type_: str) -> str:
    if type_ == "intact" or type_ == "indels":
        return type_.capitalize()
    elif type_ == "sv":
        return "SV"
    return type_


def format_result_info(path_result: Path) -> list[dict[str, str]]:
    sample_name = path_result.stem

    key_filter = ["QNAME", "NAME", "READNUM", "PERCENT"]
    key_mapping = {"QNAME": "Read ID", "READNUM": "Number of reads", "PERCENT": "Percent of reads"}
    key_order = ["Sample", "Read ID", "Label", "Allele", "Type", "Number of reads", "Percent of reads"]

    result_format = []
    for reads in io.read_jsonl(path_result):
        # Filter keys
        reads = {k: reads[k] for k in key_filter}
        # Add Label, Allele, Type from NAME
        label, allele, type_, *_ = reads["NAME"].split("_")
        reads.update({"Label": label.capitalize(), "Allele": rename_allele(allele), "Type": rename_type(type_)})
        del reads["NAME"]
        # Add Sample
        reads["Sample"] = sample_name
        # Rename keys
        for old_key, new_key in key_mapping.items():
            reads[new_key] = reads.pop(old_key)
        # Reorder keys
        reads_order = {k: reads[k] for k in key_order}
        result_format.append(reads_order)

    return result_format


def extract_all_info(batch_directory: Path) -> list[dict[str, str]]:
    results_all = []
    for path_result in sorted(batch_directory.iterdir()):
        result_format = format_result_info(path_result)
        results_all.extend(result_format)
    return results_all


def summarize_info(results_all: list[dict[str, str]]) -> list[dict[str, str]]:
    # Delete the ID column and extract unique rows
    unique_rows = []
    seen = set()  # A set to record the rows that have already been seen

    for row in results_all:
        # Delete the Read ID key
        row.pop("Read ID", None)
        # Convert to a tuple format to be able to add to the set, and add only unique rows to the seen set
        row_tuple = tuple(row.items())
        if row_tuple not in seen:
            seen.add(row_tuple)
            unique_rows.append(row)

    return sorted(unique_rows, key=lambda x: [x["Sample"], x["Label"]])


def add_allele_type(results_summary: list[dict[str, str]]) -> list[dict[str, str]]:
    for row in results_summary:
        row["Allele type"] = f"{row['Allele']} {row['Type']}"
    return results_summary


def order_allele_type(results_summary: list[dict[str, str]]) -> list[str]:
    alleles = {a["Allele"] for a in results_summary}
    alleles_insertion = {a for a in alleles if a.startswith("Insertion")}
    alleles_order = ["Control"] + sorted(alleles - alleles_insertion - {"Control"}) + sorted(alleles_insertion)

    allele_type_order = []
    for allele in alleles_order:
        for type_ in ["Intact", "Indels", "SV"]:
            allele_type_order.append(f"{allele} {type_}")

    allele_type = {a["Allele type"] for a in results_summary}
    return [at for at in allele_type_order if at in allele_type]


def output_plot(results_summary: list[dict[str, str]], report_directory: Path):
    results_plot = add_allele_type(results_summary)
    alleletype_order = order_allele_type(results_summary)

    fig = px.bar(
        results_plot,
        x="Sample",
        y="Percent of reads",
        color="Allele type",
        text="Percent of reads",
        labels={"Sample": "Samples", "Percent of reads": "Percent of reads", "Allele type": "Alelle type"},
        category_orders={"Allele type": alleletype_order},
    )
    fig.update_traces(textposition="inside", cliponaxis=False)
    fig.update_xaxes(categoryorder="category ascending")

    output_filename = Path(report_directory, "read_plot")
    fig.write_html(f"{output_filename}.html")
    # if kaleido is installed, output a pdf
    try:
        fig.write_image(f"{output_filename}.pdf")
    except Exception:
        pass


###################################################################
# main
###################################################################


def report(NAME: str) -> None:
    report_directory = Path(DAJIN_RESULTS_DIR, NAME)
    report_directory.mkdir(exist_ok=True, parents=True)
    shutil.copytree(Path(TEMP_ROOT_DIR, NAME, "report"), report_directory, dirs_exist_ok=True)

    results_all = extract_all_info(Path(TEMP_ROOT_DIR, NAME, "result"))
    results_summary = summarize_info(results_all)

    # Write to Excel
    io.write_xlsx(results_summary, Path(report_directory, "read_summary.xlsx"))

    # Write to plot as HTML and PDF
    output_plot(results_summary, report_directory)
