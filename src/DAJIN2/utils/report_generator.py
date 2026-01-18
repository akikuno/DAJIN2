from __future__ import annotations

import json
import shutil
from pathlib import Path

import plotly.express as px
import pysam

from DAJIN2.utils import fileio
from DAJIN2.utils.config import DAJIN_RESULTS_DIR, TEMP_ROOT_DIR


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
    launcher_names = ["launch_report_server.bat", "launch_report_server.command"]
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


def write_allele_viewer(report_directory: Path, genome_coordinates: dict | None) -> None:
    genome_info = None
    if genome_coordinates:
        genome = genome_coordinates.get("genome")
        chrom = genome_coordinates.get("chrom")
        start = genome_coordinates.get("start")
        end = genome_coordinates.get("end")
        if genome and chrom and start is not None and end is not None:
            genome_info = {"genome": genome, "locus": f"{chrom}:{start}-{end}"}
    genome_info_json = json.dumps(genome_info)
    viewer_template = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="utf-8">
    <title>Allele report</title>
    <script src="__IGV_URL__"></script>
    <style>
        html, body {{
            margin: 0;
            padding: 0;
            height: 100%;
            font-family: Arial, sans-serif;
            background: #f7f7f7;
        }}
        .viewer {{
            display: flex;
            flex-direction: column;
            height: 100%;
        }}
        .viewer__frame {{
            flex: 1 1 60%;
            width: 100%;
            border: none;
            background: #ffffff;
            min-height: 180px;
        }}
        .viewer__divider {{
            height: 16px;
            min-height: 16px;
            flex: 0 0 16px;
            background: #e2e8f0;
            cursor: row-resize;
            border-top: 1px solid #cbd5e0;
            border-bottom: 1px solid #cbd5e0;
            touch-action: none;
            user-select: none;
        }}
        .viewer__igv {{
            flex: 1 1 40%;
            display: flex;
            flex-direction: column;
            padding: 0.75rem 1rem 1rem;
            gap: 0.5rem;
            background: #ffffff;
            min-height: 180px;
        }}
        .viewer__title {{
            font-size: 0.95rem;
            font-weight: 600;
            color: #1a202c;
        }}
        .viewer__status {{
            font-size: 0.85rem;
            color: #4a5568;
        }}
        .viewer__igv-view {{
            flex: 1 1 auto;
            width: 100%;
            border: 1px solid #e2e8f0;
            border-radius: 0.5rem;
            overflow: hidden;
            background: #f7fafc;
        }}
    </style>
</head>
<body>
    <div class="viewer" id="viewer">
        <iframe class="viewer__frame" id="viewer-frame" title="Allele report"></iframe>
        <div class="viewer__divider" id="viewer-divider" role="separator" aria-orientation="horizontal"></div>
        <div class="viewer__igv" id="viewer-igv">
            <div class="viewer__title" id="viewer-title">Genome browser (IGV)</div>
            <div class="viewer__status" id="viewer-status">Loading...</div>
            <div class="viewer__igv-view" id="viewer-igv-view"></div>
        </div>
    </div>
    <script>
    (function() {{
        const params = new URLSearchParams(window.location.search);
        const htmlPath = params.get("html") || "";
        const bamPath = params.get("bam") || "";
        const baiPath = params.get("bai") || "";
        const vcfPath = params.get("vcf") || "";
        const sample = params.get("sample") || "";
        const title = params.get("title") || "Allele report";
        const frame = document.getElementById("viewer-frame");
        const status = document.getElementById("viewer-status");
        const igvView = document.getElementById("viewer-igv-view");
        const divider = document.getElementById("viewer-divider");
        const igvContainer = document.getElementById("viewer-igv");
        const viewer = document.getElementById("viewer");
        let splitRatio = 0.6;

        const encodePath = (path) => {
            if (!path) {
                return "";
            }
            return path
                .split("/")
                .map((segment) => encodeURIComponent(segment))
                .join("/");
        };

        if (frame) {{
            frame.src = encodePath(htmlPath);
        }}
        if (document.getElementById("viewer-title")) {{
            document.getElementById("viewer-title").textContent = title;
        }}

        const applySplit = (topHeight) => {{
            if (!viewer || !frame || !igvContainer || !divider) {{
                return;
            }}
            const total = viewer.clientHeight - divider.offsetHeight;
            const minTop = 180;
            const minBottom = 180;
            const clampedTop = Math.max(minTop, Math.min(total - minBottom, topHeight));
            frame.style.flex = "0 0 auto";
            igvContainer.style.flex = "0 0 auto";
            frame.style.height = `${{clampedTop}}px`;
            igvContainer.style.height = `${{total - clampedTop}}px`;
            splitRatio = clampedTop / total;
        }};

        if (divider && viewer) {{
            let startY = 0;
            let startTop = 0;
            let isDragging = false;
            const onMove = (event) => {{
                if (!isDragging) {{
                    return;
                }}
                event.preventDefault();
                applySplit(startTop + (event.clientY - startY));
            }};
            const stopDragging = (event) => {{
                if (!isDragging) {{
                    return;
                }}
                isDragging = false;
                document.body.style.cursor = "";
                if (event && divider.hasPointerCapture?.(event.pointerId)) {{
                    divider.releasePointerCapture(event.pointerId);
                }}
            }};
            divider.addEventListener("pointerdown", (event) => {{
                isDragging = true;
                startY = event.clientY;
                startTop = frame ? frame.getBoundingClientRect().height : 0;
                document.body.style.cursor = "row-resize";
                event.preventDefault();
                if (divider.setPointerCapture) {{
                    divider.setPointerCapture(event.pointerId);
                }}
            }});
            divider.addEventListener("pointermove", onMove);
            divider.addEventListener("pointerup", stopDragging);
            divider.addEventListener("pointercancel", stopDragging);
            window.addEventListener("resize", () => {{
                if (viewer.clientHeight > 0) {{
                    applySplit(viewer.clientHeight * splitRatio);
                }}
            }});
            if (viewer.clientHeight > 0) {{
                requestAnimationFrame(() => {{
                    applySplit(viewer.clientHeight * splitRatio);
                }});
            }}
        }}

        const hasBam = Boolean(bamPath && baiPath);
        const hasVcf = Boolean(vcfPath);
        if (!hasBam && !hasVcf) {{
            if (status) {{
                status.textContent = "No IGV tracks available.";
            }}
            return;
        }}
        if (!window.igv) {{
            if (status) {{
                status.textContent = "IGV library is not available.";
            }}
            return;
        }}

        const igvGenomeInfo = __GENOME_INFO__;
        const baseOptions = igvGenomeInfo
            ? {{ genome: igvGenomeInfo.genome, locus: igvGenomeInfo.locus }}
            : {{
                  reference: {{
                      fastaURL: sample ? encodePath(`FASTA/${{sample}}/control.fasta`) : "",
                      indexURL: sample ? encodePath(`FASTA/${{sample}}/control.fasta.fai`) : "",
                  }},
              }};
        const options = {{ ...baseOptions, tracks: [] }};
        const alignmentTrack = hasBam
            ? {{
                  name: title,
                  url: encodePath(bamPath),
                  indexURL: encodePath(baiPath),
                  indexurl: encodePath(baiPath),
                  type: "alignment",
                  format: "bam",
                  autoHeight: true,
                  viewAsPairs: true,
                  samplingDepth: 30,
                  showInsertionText: true,
                  showDeletionText: true,
              }}
            : null;
        const variantTrack = hasVcf
            ? {{
                  name: `${{title}} variants`,
                  url: encodePath(vcfPath),
                  type: "variant",
                  format: "vcf",
                  displayMode: "EXPANDED",
                  showAllSites: true,
              }}
            : null;

        igv.createBrowser(igvView, options)
            .then((browser) => {{
                const waitForRefseq = (target, timeoutMs = 2000) =>
                    new Promise((resolve) => {{
                        const start = Date.now();
                        const check = () => {{
                            const hasRefseq = Array.isArray(target?.tracks)
                                ? target.tracks.some((t) =>
                                      String(t?.name ?? "").toLowerCase().includes("refseq curated")
                                  )
                                : false;
                            if (hasRefseq) {{
                                resolve(true);
                                return;
                            }}
                            if (Date.now() - start >= timeoutMs) {{
                                resolve(false);
                                return;
                            }}
                            setTimeout(check, 100);
                        }};
                        check();
                    }});
                const enforceTrackOrder = async (target, trackConfig) => {{
                    if (!target || !Array.isArray(target.tracks)) {{
                        return;
                    }}
                    const refseqIndex = target.tracks.findIndex((t) =>
                        String(t?.name ?? "").toLowerCase().includes("refseq curated")
                    );
                    const alignmentTrackLoaded = target.tracks.find(
                        (t) => t?.type === "alignment" && String(t?.name ?? "") === String(trackConfig?.name ?? "")
                    );
                    if (refseqIndex === -1 || !alignmentTrackLoaded) {{
                        return;
                    }}
                    const alignmentIndex = target.tracks.indexOf(alignmentTrackLoaded);
                    if (alignmentIndex > refseqIndex) {{
                        return;
                    }}
                    if (typeof target.removeTrack === "function") {{
                        target.removeTrack(alignmentTrackLoaded);
                        if (typeof target.loadTrack === "function") {{
                            await target.loadTrack(trackConfig);
                        }}
                    }}
                }};
                const loadTrackSafe = async (target, trackConfig) => {{
                    if (!trackConfig || typeof target?.loadTrack !== "function") {{
                        return;
                    }}
                    await target.loadTrack(trackConfig);
                }};
                if (browser && typeof browser.loadTrack === "function") {{
                    return waitForRefseq(browser)
                        .then(() => loadTrackSafe(browser, alignmentTrack))
                        .then(() => (alignmentTrack ? enforceTrackOrder(browser, alignmentTrack) : null))
                        .then(() => loadTrackSafe(browser, variantTrack))
                        .then(() => browser);
                }}
                return browser;
            }})
            .then((browser) => {{
                if (status) {{
                    status.textContent = "";
                }}
            }})
            .catch((error) => {{
                console.error("Failed to load IGV", error);
                if (status) {{
                    status.textContent = "Failed to load IGV track.";
                }}
            }});
    }})();
    </script>
</body>
</html>
"""
    viewer_html = viewer_template.replace("{{", "{").replace("}}", "}")
    viewer_html = viewer_html.replace("__GENOME_INFO__", genome_info_json)
    viewer_html = viewer_html.replace("__IGV_URL__", "https://cdn.jsdelivr.net/npm/igv@3.7.3/dist/igv.min.js")
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


def output_plot(results_summary: list[dict[str, str]], report_directory: Path, genome_coordinates: dict | None = None):
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
    report_html_path = Path(report_directory, "report.html")
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

    base_css = """
        body {
            font-family: Arial, sans-serif;
            margin: 0;
            padding: 1.5rem;
            background-color: #f7f7f7;
        }
        #read_plot_fig {
            height: 800px !important;
        }
        .controls {
            display: flex;
            flex-wrap: wrap;
            gap: 1rem;
            align-items: flex-end;
            background-color: #ffffff;
            border-radius: 0.5rem;
            padding: 1rem;
            box-shadow: 0 1px 3px rgba(0, 0, 0, 0.1);
            margin-bottom: 1rem;
        }
        .controls label {
            display: flex;
            flex-direction: column;
            font-size: 0.9rem;
            font-weight: 600;
            gap: 0.3rem;
        }
        .controls input[type="number"] {
            padding: 0.4rem 0.6rem;
            border: 1px solid #ccc;
            border-radius: 0.3rem;
            min-width: 10rem;
        }
        .controls__actions {
            display: flex;
            gap: 0.5rem;
        }
        .controls button {
            padding: 0.5rem 0.9rem;
            border: none;
            border-radius: 0.3rem;
            font-size: 0.95rem;
            cursor: pointer;
        }
        .controls button.apply-btn {
            background: #1f77b4;
            color: #fff;
        }
        .controls button.reset-btn {
            background: #e0e0e0;
            color: #333;
        }
        .color-controls {
            background-color: #ffffff;
            border-radius: 0.5rem;
            padding: 1rem;
            box-shadow: 0 1px 3px rgba(0, 0, 0, 0.1);
            margin-bottom: 1.5rem;
        }
        .color-controls__title {
            font-size: 0.95rem;
            font-weight: 600;
            margin-bottom: 0.75rem;
        }
        .color-controls__list {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(220px, 1fr));
            gap: 0.75rem;
        }
        .export-controls {
            display: flex;
            flex-wrap: wrap;
            gap: 0.5rem;
        }
        .export-controls button {
            padding: 0.5rem 0.9rem;
            border: none;
            border-radius: 0.3rem;
            font-size: 0.95rem;
            cursor: pointer;
            background: #4a5568;
            color: #fff;
        }
        .color-picker-item {
            display: flex;
            justify-content: space-between;
            align-items: center;
            gap: 0.75rem;
            padding: 0.5rem 0.75rem;
            border: 1px solid #e0e0e0;
            border-radius: 0.5rem;
            background: #fafafa;
        }
        .color-picker-item span {
            font-size: 0.95rem;
            font-weight: 600;
            cursor: pointer;
        }
        .color-picker-item input[type="color"] {
            width: 2.4rem;
            height: 2.4rem;
            padding: 0;
            border: none;
            background: transparent;
            cursor: pointer;
        }
        .plot-container {
            background-color: #ffffff;
            border-radius: 0.5rem;
            padding: 1rem;
            box-shadow: 0 1px 3px rgba(0, 0, 0, 0.1);
        }
        @media (max-width: 600px) {
            body {
                padding: 1rem;
            }
            .controls {
                flex-direction: column;
                align-items: stretch;
            }
            .controls__actions {
                justify-content: flex-start;
            }
        }
    """
    controls_block = f"""
    <div class="controls">
        <label>
            Font size
            <input id="plot-font-size" type="number" min="6" max="72" value="{int(initial_font_size)}">
        </label>
        <div class="controls__actions">
            <button class="apply-btn" id="plot-style-apply">Apply</button>
            <button class="reset-btn" id="plot-style-reset">Reset</button>
        </div>
    </div>
    <div class="color-controls">
        <div class="color-controls__title">Bar colors (click label to open picker)</div>
        <div class="color-controls__list" id="color-picker-list"></div>
    </div>
    <div class="export-controls" id="export-controls"></div>
"""

    script_block = f"""
    <script>
    (function() {{
        const figure = document.getElementById("{div_id}");
        if (!figure) {{
            return;
        }}

        const fontInput = document.getElementById("plot-font-size");
        const applyButton = document.getElementById("plot-style-apply");
        const resetButton = document.getElementById("plot-style-reset");
        const colorPickerList = document.getElementById("color-picker-list");
        const exportControls = document.getElementById("export-controls");

        const clone = (obj) => JSON.parse(JSON.stringify(obj ?? (Array.isArray(obj) ? [] : {{}})));
        const initialData = clone(figure.data);
        const initialLayout = clone(figure.layout);
        const initialFontValue = fontInput.value;

        const traceGroups = {trace_groups_json};
        const exportButtons = {export_buttons_json};

        const ensureHex = (color) => {{
            if (typeof color !== "string") {{
                return "#1f77b4";
            }}
            if (color.startsWith("#")) {{
                if (color.length === 4) {{
                    return "#" + color[1] + color[1] + color[2] + color[2] + color[3] + color[3];
                }}
                return color;
            }}
            const rgbMatch = color.match(/rgba?\\(([^)]+)\\)/i);
            if (rgbMatch) {{
                const parts = rgbMatch[1].split(",").map((item) => parseFloat(item.trim()));
                const toHex = (value) => {{
                    const v = Math.max(0, Math.min(255, Math.round(value)));
                    return v.toString(16).padStart(2, "0");
                }};
                return "#" + toHex(parts[0]) + toHex(parts[1]) + toHex(parts[2]);
            }}
            return "#1f77b4";
        }};

        const getTraceColor = (trace) => {{
            if (!trace || !trace.marker) {{
                return "#1f77b4";
            }}
            const color = trace.marker.color;
            if (Array.isArray(color) && color.length > 0) {{
                return ensureHex(color[0]);
            }}
            if (typeof color === "string") {{
                return ensureHex(color);
            }}
            return "#1f77b4";
        }};

        const buildColorPickers = () => {{
            colorPickerList.innerHTML = "";
            Object.entries(traceGroups).forEach(([name, indices]) => {{
                const trace = figure.data[indices[0]];

                const wrapper = document.createElement("label");
                wrapper.className = "color-picker-item";

                const title = document.createElement("span");
                title.textContent = name;

                const picker = document.createElement("input");
                picker.type = "color";
                picker.value = getTraceColor(trace);

                picker.addEventListener("input", () => {{
                    Plotly.restyle(figure, {{"marker.color": picker.value}}, indices);
                }});
                title.addEventListener("click", () => picker.click());

                wrapper.appendChild(title);
                wrapper.appendChild(picker);
                colorPickerList.appendChild(wrapper);
            }});
        }};

        const applyFontSize = () => {{
            const fontSize = parseInt(fontInput.value, 10);
            if (Number.isNaN(fontSize) || fontSize <= 0) {{
                return;
            }}
            Plotly.relayout(figure, {{
                "font.size": fontSize,
                "xaxis.title.font.size": fontSize,
                "xaxis.tickfont.size": fontSize,
                "yaxis.title.font.size": fontSize,
                "yaxis.tickfont.size": fontSize,
                "legend.font.size": fontSize
            }});
            Plotly.restyle(figure, {{"textfont.size": fontSize}});
        }};

        const resetStyle = () => {{
            Plotly.react(figure, clone(initialData), clone(initialLayout));
            fontInput.value = initialFontValue;
            buildColorPickers();
        }};

        const buildExportButtons = () => {{
            exportControls.innerHTML = "";
            exportButtons.forEach((btn) => {{
                const button = document.createElement("button");
                button.textContent = btn.label;
                button.dataset.format = btn.format;
                button.addEventListener("click", async () => {{
                    try {{
                        await Plotly.downloadImage(figure, {{
                            format: btn.format,
                            filename: "read_plot",
                            scale: 600 / 96,
                        }});
                    }} catch (error) {{
                        console.error("Failed to export image", error);
                        alert("Failed to export " + btn.format.toUpperCase() + " file. Please ensure Plotly is fully loaded.");
                    }}
                }});
                exportControls.appendChild(button);
            }});
        }};

        applyButton.addEventListener("click", applyFontSize);
        resetButton.addEventListener("click", resetStyle);
        buildColorPickers();
        buildExportButtons();
    }})();
    </script>
"""

    report_hint_block = """
    <div class="report-hint">
        Click an allele segment to open the detailed report and genome browser.
    </div>
    """

    modal_block = """
    <div class="modal" id="allele-modal" aria-hidden="true">
        <div class="modal__backdrop" data-modal-close></div>
        <div class="modal__content" role="dialog" aria-modal="true" aria-labelledby="allele-modal-title">
            <div class="modal__header">
                <div class="modal__title" id="allele-modal-title">Allele report</div>
                <button class="modal__close" id="allele-modal-close" type="button">Close</button>
            </div>
            <div class="modal__body">
                <iframe class="modal__frame" id="allele-modal-frame" title="Allele report"></iframe>
                <div class="modal__divider" id="allele-modal-divider" role="separator" aria-orientation="horizontal"></div>
                <div class="modal__igv">
                    <div class="modal__igv-title">Genome browser (IGV)</div>
                    <div class="modal__igv-status" id="allele-igv-status">Waiting for allele selection.</div>
                    <div class="modal__igv-view" id="allele-igv-view"></div>
                </div>
            </div>
            <div class="modal__footer">
                <a class="modal__link" id="allele-modal-link" target="_blank" rel="noopener">Open in new tab</a>
            </div>
        </div>
    </div>
    """

    modal_css = """
        .report-hint {
            margin-bottom: 1rem;
            font-size: 0.95rem;
            color: #2d3748;
        }
        .modal {
            position: fixed;
            inset: 0;
            display: none;
            align-items: center;
            justify-content: center;
            padding: 1.5rem;
            z-index: 1000;
        }
        .modal.is-open {
            display: flex;
        }
        .modal__backdrop {
            position: absolute;
            inset: 0;
            background: rgba(0, 0, 0, 0.45);
        }
        .modal__content {
            position: relative;
            background: #ffffff;
            border-radius: 0.75rem;
            width: min(1100px, 94vw);
            height: 90vh;
            max-height: 90vh;
            display: flex;
            flex-direction: column;
            box-shadow: 0 12px 30px rgba(0, 0, 0, 0.2);
            overflow: hidden;
            z-index: 1;
        }
        .modal__body {
            display: flex;
            flex-direction: column;
            background: #ffffff;
            flex: 1 1 auto;
            min-height: 0;
            overflow: hidden;
        }
        .modal__header,
        .modal__footer {
            display: flex;
            justify-content: space-between;
            align-items: center;
            gap: 1rem;
            padding: 0.9rem 1.2rem;
            background: #f7fafc;
            border-bottom: 1px solid #e2e8f0;
        }
        .modal__footer {
            border-top: 1px solid #e2e8f0;
            border-bottom: none;
        }
        .modal__title {
            font-size: 1rem;
            font-weight: 600;
            color: #1a202c;
        }
        .modal__close {
            border: none;
            background: #e2e8f0;
            color: #1a202c;
            padding: 0.4rem 0.8rem;
            border-radius: 0.4rem;
            cursor: pointer;
        }
        .modal__frame {
            width: 100%;
            height: 100%;
            border: none;
            background: #ffffff;
            flex: 1 1 60%;
            min-height: 180px;
        }
        .modal__divider {
            height: 16px;
            min-height: 16px;
            flex: 0 0 16px;
            background: #e2e8f0;
            cursor: row-resize;
            border-top: 1px solid #cbd5e0;
            border-bottom: 1px solid #cbd5e0;
            touch-action: none;
            user-select: none;
        }
        .modal__igv {
            display: flex;
            flex-direction: column;
            gap: 0.5rem;
            padding: 0.75rem 1rem 1rem;
            border-top: 1px solid #e2e8f0;
            background: #ffffff;
            flex: 1 1 40%;
            min-height: 180px;
        }
        .modal__igv-title {
            font-size: 0.9rem;
            font-weight: 600;
            color: #1a202c;
        }
        .modal__igv-status {
            font-size: 0.85rem;
            color: #4a5568;
        }
        .modal__igv-view {
            width: 100%;
            height: 100%;
            border: 1px solid #e2e8f0;
            border-radius: 0.5rem;
            overflow: hidden;
            background: #f7fafc;
        }
        .modal__link {
            text-decoration: none;
            color: #3182ce;
            font-weight: 600;
        }
    """

    report_style_block = build_style_block(base_css, modal_css)

    genome_info = None
    if genome_coordinates:
        genome = genome_coordinates.get("genome")
        chrom = genome_coordinates.get("chrom")
        start = genome_coordinates.get("start")
        end = genome_coordinates.get("end")
        if genome and chrom and start is not None and end is not None:
            genome_info = {"genome": genome, "locus": f"{chrom}:{start}-{end}"}
    genome_info_json = json.dumps(genome_info)

    igv_lib_block = """
    <script src="https://cdn.jsdelivr.net/npm/igv@3.7.3/dist/igv.min.js"></script>
    """

    report_script_block = """
    <script>
    (function() {
        const figure = document.getElementById("read_plot_fig");
        const modal = document.getElementById("allele-modal");
        const modalTitle = document.getElementById("allele-modal-title");
        const modalFrame = document.getElementById("allele-modal-frame");
        const modalLink = document.getElementById("allele-modal-link");
        const modalClose = document.getElementById("allele-modal-close");
        const modalBackdrop = modal ? modal.querySelector("[data-modal-close]") : null;
        const modalBody = modal ? modal.querySelector(".modal__body") : null;
        const modalDivider = document.getElementById("allele-modal-divider");
        const igvContainer = modal ? modal.querySelector(".modal__igv") : null;
        const igvStatus = document.getElementById("allele-igv-status");
        let igvView = document.getElementById("allele-igv-view");
        const igvGenomeInfo = __GENOME_INFO__;
        let igvBrowser = null;
        let igvRequestId = 0;
        let splitRatio = 0.6;

        if (!figure || !modal) {
            return;
        }

        const formatPercent = (value) => {
            if (value === null || value === undefined || value === "") {
                return "";
            }
            return `${value}%`;
        };

        const encodePath = (path) => {
            if (!path) {
                return "";
            }
            return path
                .split("/")
                .map((segment) => encodeURIComponent(segment))
                .join("/");
        };

        const formatAlleleType = (allele, type) => {
            const alleleText = String(allele ?? "").trim();
            const typeText = String(type ?? "").trim();
            const typeLower = typeText.toLowerCase();
            if (!typeLower || typeLower === "intact") {
                return alleleText;
            }
            if (typeLower === "indels") {
                return `${alleleText} with indels`.trim();
            }
            return `${alleleText} ${typeText}`.trim();
        };

        const buildTitle = (point) => {
            const custom = Array.isArray(point.customdata) ? point.customdata : [];
            const label = custom[1];
            const allele = custom[2];
            const type = custom[3];
            const percent = formatPercent(custom[4]);
            const alleleType = formatAlleleType(allele, type);
            const labelText = label ? label.toLowerCase() : "";
            const percentText = percent ? `(${percent})` : "";
            const parts = [
                point.x,
                labelText ? `${labelText}` : "",
                alleleType || "",
                percentText,
            ].filter(Boolean);
            return parts.join(" ").replace(/\s+/g, " ").trim();
        };

        const buildIgvPaths = (path) => {
            if (!path) {
                return null;
            }
            const parts = path.split("/");
            if (parts.length < 3) {
                return null;
            }
            const htmlIndex = parts.indexOf("HTML");
            const baseIndex = htmlIndex >= 0 ? htmlIndex : 0;
            if (parts.length <= baseIndex + 2) {
                return null;
            }
            const sample = parts[baseIndex + 1];
            const filename = parts[parts.length - 1];
            const stem = filename.replace(/\\.html?$/i, "");
            const prefix = `${sample}_`;
            const header = stem.startsWith(prefix) ? stem.slice(prefix.length) : stem;
            return {
                sample,
                bam: [".igvjs", sample, `${header}.bam`].join("/"),
                bai: [".igvjs", sample, `${header}.bam.bai`].join("/"),
                vcf: ["VCF", sample, `${sample}_${header}.vcf`].join("/"),
            };
        };

        const buildViewerUrl = (path, title) => {
            const igvPaths = buildIgvPaths(path);
            if (!igvPaths) {
                return encodePath(path);
            }
            const params = new URLSearchParams({
                html: path,
                bam: igvPaths.bam,
                bai: igvPaths.bai,
                vcf: igvPaths.vcf,
                sample: igvPaths.sample,
                title: title || "",
            });
            return `allele_viewer.html?${params.toString()}`;
        };

        const buildIgvBaseOptions = (sample) => {
            const baseOptions = igvGenomeInfo
                ? { genome: igvGenomeInfo.genome, locus: igvGenomeInfo.locus }
                : {
                      reference: {
                          fastaURL: encodePath(`FASTA/${sample}/control.fasta`),
                          indexURL: encodePath(`FASTA/${sample}/control.fasta.fai`),
                      },
                  };
            return { ...baseOptions, tracks: [] };
        };

        const buildIgvTrack = (bamUrl, baiUrl, trackName) => ({
            name: trackName || "Allele",
            url: bamUrl,
            indexURL: baiUrl,
            indexurl: baiUrl,
            type: "alignment",
            format: "bam",
            autoHeight: true,
            viewAsPairs: true,
                        samplingDepth: 30,
            showInsertionText: true,
            showDeletionText: true,
        });

        const buildVcfTrack = (vcfUrl, trackName) => ({
            name: trackName ? `${trackName} variants` : "Variants",
            url: vcfUrl,
            type: "variant",
            format: "vcf",
            displayMode: "EXPANDED",
            showAllSites: true,
        });

        const resetIgvView = () => {
            if (!igvContainer) {
                return;
            }
            if (igvView && igvView.parentNode) {
                igvView.parentNode.removeChild(igvView);
            }
            igvView = document.createElement("div");
            igvView.id = "allele-igv-view";
            igvView.className = "modal__igv-view";
            igvContainer.appendChild(igvView);
        };

        const updateIgv = async (path, title) => {
            if (!igvView || !igvStatus) {
                return;
            }
            const igvPaths = buildIgvPaths(path);
            if (!igvPaths) {
                igvStatus.textContent = "No IGV tracks available for this allele.";
                igvView.innerHTML = "";
                return;
            }
            if (!window.igv) {
                igvStatus.textContent = "IGV library is not available.";
                igvView.innerHTML = "";
                return;
            }
            const hasBam = Boolean(igvPaths.bam && igvPaths.bai);
            const hasVcf = Boolean(igvPaths.vcf);
            if (!hasBam && !hasVcf) {
                igvStatus.textContent = "No IGV tracks available for this allele.";
                igvView.innerHTML = "";
                return;
            }
            igvStatus.textContent = "Loading IGV...";
            const options = buildIgvBaseOptions(igvPaths.sample);
            const currentRequestId = ++igvRequestId;
            try {
                const addCacheBuster = (url, seed) => {
                    if (!url) {
                        return url;
                    }
                    const joiner = url.includes("?") ? "&" : "?";
                    return `${url}${joiner}v=${seed}`;
                };
                const track = hasBam ? buildIgvTrack(igvPaths.bam, igvPaths.bai, title) : null;
                if (track) {
                    const encodedTrackUrl = addCacheBuster(encodePath(track.url), currentRequestId);
                    const encodedIndexUrl = addCacheBuster(encodePath(track.indexURL), currentRequestId);
                    track.url = encodedTrackUrl;
                    track.indexURL = encodedIndexUrl;
                    track.indexurl = encodedIndexUrl;
                }
                const variantTrack = hasVcf ? buildVcfTrack(igvPaths.vcf, title) : null;
                if (variantTrack) {
                    variantTrack.url = addCacheBuster(encodePath(variantTrack.url), currentRequestId);
                }
                const nextSignature = [igvPaths.bam, igvPaths.bai, igvPaths.vcf, title || ""].join("|");
                const previousSignature = igvView.dataset.igvSignature || "";
                igvView.dataset.igvSignature = nextSignature;
                if (igvBrowser && typeof igvBrowser.dispose === "function") {
                    igvBrowser.dispose();
                } else if (igvBrowser && window.igv && typeof igv.removeBrowser === "function") {
                    igv.removeBrowser(igvBrowser);
                }
                igvBrowser = null;
                resetIgvView();
                const browser = await igv.createBrowser(igvView, options);
                if (currentRequestId !== igvRequestId) {
                    if (browser && typeof browser.dispose === "function") {
                        browser.dispose();
                    }
                    return;
                }
                const waitForRefseq = (target, timeoutMs = 2000) =>
                    new Promise((resolve) => {
                        const start = Date.now();
                        const check = () => {
                            const hasRefseq = Array.isArray(target?.tracks)
                                ? target.tracks.some((t) =>
                                      String(t?.name ?? "").toLowerCase().includes("refseq curated")
                                  )
                                : false;
                            if (hasRefseq) {
                                resolve(true);
                                return;
                            }
                            if (Date.now() - start >= timeoutMs) {
                                resolve(false);
                                return;
                            }
                            setTimeout(check, 100);
                        };
                        check();
                    });
                if (browser && typeof browser.loadTrack === "function") {
                    await waitForRefseq(browser);
                    if (track) {
                        await browser.loadTrack(track);
                        await enforceTrackOrder(browser, track);
                    }
                    if (variantTrack) {
                        await browser.loadTrack(variantTrack);
                    }
                }
                igvBrowser = browser;
                window.__igvBrowser = browser;
                window.__igvLastOptions = options;
                igvStatus.textContent = "";
            } catch (error) {
                console.error("Failed to load IGV", error);
                igvStatus.textContent = "Failed to load IGV track.";
            }
        };

        const enforceTrackOrder = async (browser, trackConfig) => {
            if (!browser || !Array.isArray(browser.tracks)) {
                return;
            }
            const refseqIndex = browser.tracks.findIndex((t) =>
                String(t?.name ?? "").toLowerCase().includes("refseq curated")
            );
            const alignmentTrack = browser.tracks.find(
                (t) => t?.type === "alignment" && String(t?.name ?? "") === String(trackConfig?.name ?? "")
            );
            if (refseqIndex === -1 || !alignmentTrack) {
                return;
            }
            const alignmentIndex = browser.tracks.indexOf(alignmentTrack);
            if (alignmentIndex > refseqIndex) {
                return;
            }
            if (typeof browser.removeTrack === "function") {
                browser.removeTrack(alignmentTrack);
                if (typeof browser.loadTrack === "function") {
                    await browser.loadTrack(trackConfig);
                }
            }
        };

        const applySplit = (topHeight) => {
            if (!modalBody || !modalFrame || !igvContainer || !modalDivider) {
                return;
            }
            const total = modalBody.clientHeight - modalDivider.offsetHeight;
            if (total <= 0) {
                return;
            }
            const minTop = 180;
            const minBottom = 180;
            const clampedTop = Math.max(minTop, Math.min(total - minBottom, topHeight));
            modalFrame.style.flex = "0 0 auto";
            igvContainer.style.flex = "0 0 auto";
            modalFrame.style.height = `${clampedTop}px`;
            igvContainer.style.height = `${total - clampedTop}px`;
            splitRatio = clampedTop / total;
        };

        const initDivider = () => {
            if (!modalDivider || !modalBody || !modalFrame) {
                return;
            }
            let startY = 0;
            let startTop = 0;
            let isDragging = false;
            const onMove = (event) => {
                if (!isDragging) {
                    return;
                }
                event.preventDefault();
                applySplit(startTop + (event.clientY - startY));
            };
            const stopDragging = (event) => {
                if (!isDragging) {
                    return;
                }
                isDragging = false;
                document.body.style.cursor = "";
                if (event && modalDivider.hasPointerCapture?.(event.pointerId)) {
                    modalDivider.releasePointerCapture(event.pointerId);
                }
            };
            modalDivider.addEventListener("pointerdown", (event) => {
                isDragging = true;
                startY = event.clientY;
                startTop = modalFrame.getBoundingClientRect().height;
                document.body.style.cursor = "row-resize";
                event.preventDefault();
                if (modalDivider.setPointerCapture) {
                    modalDivider.setPointerCapture(event.pointerId);
                }
            });
            modalDivider.addEventListener("pointermove", onMove);
            modalDivider.addEventListener("pointerup", stopDragging);
            modalDivider.addEventListener("pointercancel", stopDragging);
            window.addEventListener("resize", () => {
                if (modalBody.clientHeight > 0) {
                    applySplit(modalBody.clientHeight * splitRatio);
                }
            });
        };

        const openModal = (path, title) => {
            modal.classList.add("is-open");
            modal.setAttribute("aria-hidden", "false");
            const encodedPath = encodePath(path);
            modalFrame.src = encodedPath;
            modalTitle.textContent = title || "Allele report";
            modalLink.href = buildViewerUrl(path, title);
            updateIgv(path, title);
            if (modalBody) {
                requestAnimationFrame(() => {
                    if (modalBody.clientHeight > 0) {
                        applySplit(modalBody.clientHeight * splitRatio);
                    }
                });
            }
        };

        const closeModal = () => {
            modal.classList.remove("is-open");
            modal.setAttribute("aria-hidden", "true");
            modalFrame.src = "";
            modalTitle.textContent = "Allele report";
            modalLink.removeAttribute("href");
            if (igvView) {
                igvView.innerHTML = "";
            }
            if (igvStatus) {
                igvStatus.textContent = "Waiting for allele selection.";
            }
            igvBrowser = null;
        };

        initDivider();

        modalClose.addEventListener("click", closeModal);
        if (modalBackdrop) {
            modalBackdrop.addEventListener("click", closeModal);
        }

        figure.on("plotly_click", (event) => {
            if (!event || !event.points || !event.points.length) {
                return;
            }
            const point = event.points[0];
            const custom = Array.isArray(point.customdata) ? point.customdata : [];
            const detailPath = custom[0];
            if (!detailPath) {
                alert("No detailed report found for this allele.");
                return;
            }
            openModal(detailPath, buildTitle(point));
        });
    })();
    </script>
    """

    report_script_block = report_script_block.replace("__GENOME_INFO__", genome_info_json)

    inject_plot_assets(
        report_html_path,
        report_style_block,
        [report_hint_block, controls_block, modal_block],
        [script_block, igv_lib_block, report_script_block],
    )


###################################################################
# main
###################################################################


def report(NAME: str) -> None:
    report_directory = Path(DAJIN_RESULTS_DIR, NAME)
    report_directory.mkdir(exist_ok=True, parents=True)
    shutil.copytree(Path(TEMP_ROOT_DIR, NAME, "report"), report_directory, dirs_exist_ok=True)
    copy_report_launchers(report_directory)

    results_all = extract_all_info(Path(TEMP_ROOT_DIR, NAME, "result"))
    results_summary = summarize_info(results_all)

    # Write to Excel
    fileio.write_xlsx(results_summary, Path(report_directory, "read_summary.xlsx"))

    # Write to plot as HTML and PDF
    genome_coordinates = load_genome_coordinates(NAME)
    if genome_coordinates is None:
        ensure_reference_indexes(report_directory)
    write_allele_viewer(report_directory, genome_coordinates)
    output_plot(results_summary, report_directory, genome_coordinates)
