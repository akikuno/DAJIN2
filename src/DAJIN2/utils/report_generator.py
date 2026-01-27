from __future__ import annotations

import json
import os
import shutil
from pathlib import Path

import plotly.express as px
import pysam

from DAJIN2.utils import fileio
from DAJIN2.utils.config import DAJIN_RESULTS_DIR, TEMP_ROOT_DIR

TEMPLATES_DIR = Path(__file__).resolve().parent.parent / "templates"


def read_template(name: str) -> str:
    return Path(TEMPLATES_DIR, name).read_text(encoding="utf-8")


def format_result_info(path_result: Path) -> list[dict[str, str]]:
    sample_name = path_result.stem

    key_filter = ["QNAME", "NAME", "READNUM", "PERCENT"]
    key_mapping = {"QNAME": "Read ID", "READNUM": "Number of reads", "PERCENT": "Percent of reads"}
    key_order = ["Sample", "Read ID", "Label", "Allele", "Type", "Number of reads", "Percent of reads"]

    result_format = []
    for reads in fileio.read_jsonl(path_result):
        # Filter keys
        reads = {k: reads[k] for k in key_filter}
        # Add Label, Allele, Type from NAME
        label, allele, type_, *_ = reads["NAME"].split("_")
        reads.update({"Label": label.capitalize(), "Allele": allele, "Type": type_})
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
        row["Allele type"] = format_allele_type_label(row.get("Allele", ""), row.get("Type", ""))
    return results_summary


def format_allele_type_label(allele: str, type_: str) -> str:
    allele_text = str(allele).strip()
    type_text = str(type_).strip()
    type_lower = type_text.lower()
    if type_lower in {"", "intact"}:
        return allele_text
    if type_lower == "indels":
        return f"{allele_text} with indels".strip()
    return f"{allele_text} {type_text}".strip()


def load_genome_coordinates(report_name: str) -> dict | None:
    path_coordinates = Path(TEMP_ROOT_DIR, report_name, "cache", "genome_coordinates.jsonl")
    if not path_coordinates.exists():
        return None
    try:
        return next(fileio.read_jsonl(path_coordinates))
    except StopIteration:
        return None


def ensure_reference_indexes(report_directory: Path) -> None:
    fasta_root = Path(report_directory, "FASTA")
    if not fasta_root.exists():
        return
    for control_fasta in fasta_root.glob("*/control.fasta"):
        path_fai = Path(str(control_fasta) + ".fai")
        if not path_fai.exists():
            pysam.faidx(str(control_fasta))


def copy_report_launchers(report_directory: Path) -> None:
    source_root = Path(__file__).resolve().parent.parent
    launcher_names = ["launch_report_windows.bat", "launch_report_mac.command"]
    for name in launcher_names:
        source_path = Path(source_root, name)
        if not source_path.exists():
            continue
        target_path = Path(report_directory, name)
        shutil.copy2(source_path, target_path)
        if target_path.suffix == ".command":
            try:
                target_path.chmod(0o755)
            except PermissionError:
                pass


def write_allele_viewer(report_directory: Path, genome_coordinates: dict | None, asset_prefix: str = "") -> None:
    genome_info = None
    if genome_coordinates:
        genome = genome_coordinates.get("genome")
        chrom = genome_coordinates.get("chrom")
        start = genome_coordinates.get("start")
        end = genome_coordinates.get("end")
        if genome and chrom and start is not None and end is not None:
            genome_info = {"genome": genome, "locus": f"{chrom}:{start}-{end}"}
    genome_info_json = json.dumps(genome_info)
    viewer_html = read_template("allele_viewer.html")
    viewer_css = read_template("allele_viewer.css")
    igv_helpers_js = read_template("igv_helpers.js")
    viewer_js = read_template("allele_viewer.js")
    viewer_js = "\n".join([igv_helpers_js, viewer_js])
    viewer_html = viewer_html.replace("__ALLELE_VIEWER_CSS__", viewer_css)
    viewer_html = viewer_html.replace("__ALLELE_VIEWER_JS__", viewer_js)
    viewer_html = viewer_html.replace("__GENOME_INFO__", genome_info_json)
    viewer_html = viewer_html.replace("__IGV_URL__", "https://cdn.jsdelivr.net/npm/igv@3.7.3/dist/igv.min.js")
    viewer_html = viewer_html.replace("__ASSET_PREFIX__", json.dumps(asset_prefix))
    Path(report_directory, "allele_viewer.html").write_text(viewer_html, encoding="utf-8")


def build_html_lookup(report_directory: Path) -> dict[tuple[str, str, str, str], str]:
    html_root = Path(report_directory, "HTML")
    if not html_root.exists():
        return {}

    lookup: dict[tuple[str, str, str, str], str] = {}
    for sample_dir in sorted(html_root.iterdir()):
        if not sample_dir.is_dir():
            continue
        sample_name = sample_dir.name
        sample_prefix = f"{sample_name}_"
        for html_file in sorted(sample_dir.glob("*.html")):
            stem = html_file.stem
            if not stem.lower().startswith(sample_prefix.lower()):
                continue
            header = stem[len(sample_prefix) :]
            parts = header.split("_")
            if len(parts) < 3:
                continue
            label = parts[0]
            allele = parts[1]
            type_ = parts[2]
            key = (sample_name.lower(), label.lower(), allele.lower(), type_.lower())
            lookup.setdefault(key, html_file.relative_to(report_directory).as_posix())

    return lookup


def attach_html_paths(results_summary: list[dict[str, str]], report_directory: Path) -> list[dict[str, str]]:
    html_lookup = build_html_lookup(report_directory)
    for row in results_summary:
        key = (
            str(row.get("Sample", "")).lower(),
            str(row.get("Label", "")).lower(),
            str(row.get("Allele", "")).lower(),
            str(row.get("Type", "")).lower(),
        )
        row["HTML path"] = html_lookup.get(key, "")
    return results_summary


def order_allele_type(results_summary: list[dict[str, str]]) -> list[str]:
    alleles = {a["Allele"] for a in results_summary}
    alleles_insertion = {a for a in alleles if a.startswith("Insertion")}
    alleles_order = ["Control"] + sorted(alleles - alleles_insertion - {"Control"}) + sorted(alleles_insertion)

    allele_type_order = []
    for allele in alleles_order:
        for type_ in ["Intact", "Indels", "SV"]:
            allele_type_order.append(format_allele_type_label(allele, type_))

    allele_type = {a["Allele type"] for a in results_summary}
    return [at for at in allele_type_order if at in allele_type]


def build_style_block(base_css: str, extra_css: str = "") -> str:
    return f"<style>\n{base_css}\n{extra_css}\n</style>"


def inject_plot_assets(
    html_path: Path,
    style_block: str,
    body_blocks: list[str],
    script_blocks: list[str],
) -> None:
    html_content = html_path.read_text(encoding="utf-8")
    if "</head>" in html_content:
        html_content = html_content.replace("</head>", f"{style_block}\n</head>", 1)
    if "<body>" in html_content:
        html_content = html_content.replace("<body>", "<body>\n", 1)
    body_injection = "\n".join(body_blocks + script_blocks)
    html_content = html_content.replace("</body>", f"{body_injection}\n</body>", 1)
    html_path.write_text(html_content, encoding="utf-8")


def output_plot(
    results_summary: list[dict[str, str]],
    report_directory: Path,
    report_html_directory: Path,
    genome_coordinates: dict | None = None,
    asset_prefix: str = "",
):
    results_plot = add_allele_type(results_summary)
    results_plot = attach_html_paths(results_plot, report_directory)
    alleletype_order = order_allele_type(results_summary)

    fig = px.bar(
        results_plot,
        x="Sample",
        y="Percent of reads",
        color="Allele type",
        text="Percent of reads",
        labels={"Sample": "Samples", "Percent of reads": "Percent of reads", "Allele type": "Alelle type"},
        category_orders={"Allele type": alleletype_order},
        custom_data=["HTML path", "Label", "Allele", "Type", "Percent of reads"],
    )
    fig.update_traces(textposition="inside", cliponaxis=False)
    fig.update_xaxes(categoryorder="category ascending")

    div_id = "read_plot_fig"
    report_html_path = Path(report_html_directory, "report.html")
    fig.write_html(report_html_path, include_plotlyjs="cdn", full_html=True, div_id=div_id)

    trace_groups: dict[str, list[int]] = {}
    for index, trace in enumerate(fig.data):
        name = trace.name or f"Series {index + 1}"
        trace_groups.setdefault(name, []).append(index)

    initial_font_size = fig.layout.font.size or 18
    trace_groups_json = json.dumps(trace_groups)
    export_buttons = [
        {"format": "png", "label": "Download PNG"},
        {"format": "svg", "label": "Download SVG"},
    ]
    export_buttons_json = json.dumps(export_buttons)

    base_css = read_template("report_base.css")
    modal_css = read_template("report_modal.css")
    report_style_block = build_style_block(base_css, modal_css)

    controls_block = read_template("report_controls.html").replace(
        "__INITIAL_FONT_SIZE__",
        str(int(initial_font_size)),
    )
    script_js = read_template("report_controls.js")
    script_js = script_js.replace("__PLOT_DIV_ID__", div_id)
    script_js = script_js.replace("__TRACE_GROUPS__", trace_groups_json)
    script_js = script_js.replace("__EXPORT_BUTTONS__", export_buttons_json)
    script_block = "<script>\n" + script_js + "\n</script>"

    modal_block = read_template("report_modal.html")

    genome_info = None
    if genome_coordinates:
        genome = genome_coordinates.get("genome")
        chrom = genome_coordinates.get("chrom")
        start = genome_coordinates.get("start")
        end = genome_coordinates.get("end")
        if genome and chrom and start is not None and end is not None:
            genome_info = {"genome": genome, "locus": f"{chrom}:{start}-{end}"}
    genome_info_json = json.dumps(genome_info)

    igv_helpers_js = read_template("igv_helpers.js")
    report_script_js = read_template("report_main.js")
    report_script_js = "\n".join([igv_helpers_js, report_script_js])
    report_script_js = report_script_js.replace("__GENOME_INFO__", genome_info_json)
    report_script_js = report_script_js.replace("__ASSET_PREFIX__", json.dumps(asset_prefix))
    report_script_block = "<script>\n" + report_script_js + "\n</script>"


    inject_plot_assets(
        report_html_path,
        report_style_block,
        [controls_block, modal_block],
        [script_block, report_script_block],
    )


###################################################################
# main
###################################################################


def report(NAME: str) -> None:
    report_directory = Path(DAJIN_RESULTS_DIR, NAME)
    report_directory.mkdir(exist_ok=True, parents=True)
    shutil.copytree(Path(TEMP_ROOT_DIR, NAME, "report"), report_directory, dirs_exist_ok=True)
    report_html_directory = Path(report_directory, ".report")
    report_html_directory.mkdir(exist_ok=True)
    copy_report_launchers(report_directory)

    results_all = extract_all_info(Path(TEMP_ROOT_DIR, NAME, "result"))
    results_summary = summarize_info(results_all)

    # Write to Excel
    fileio.write_xlsx(results_summary, Path(report_directory, "read_summary.xlsx"))

    # Write to plot as HTML and PDF
    genome_coordinates = load_genome_coordinates(NAME)
    if genome_coordinates is None:
        ensure_reference_indexes(report_directory)
    asset_prefix = Path(os.path.relpath(report_directory, report_html_directory)).as_posix()
    if asset_prefix == ".":
        asset_prefix = ""
    write_allele_viewer(report_html_directory, genome_coordinates, asset_prefix)
    output_plot(results_summary, report_directory, report_html_directory, genome_coordinates, asset_prefix)
