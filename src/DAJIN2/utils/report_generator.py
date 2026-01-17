from __future__ import annotations

import json
import shutil
from pathlib import Path

import plotly.express as px

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
        row["Allele type"] = f"{row['Allele']} {row['Type']}"
    return results_summary


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
            allele_type_order.append(f"{allele} {type_}")

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


def output_plot(results_summary: list[dict[str, str]], report_directory: Path):
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
        Click an allele segment to open the detailed report.
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
            <iframe class="modal__frame" id="allele-modal-frame" title="Allele report"></iframe>
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
            max-height: 90vh;
            display: flex;
            flex-direction: column;
            box-shadow: 0 12px 30px rgba(0, 0, 0, 0.2);
            overflow: hidden;
            z-index: 1;
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
            height: 70vh;
            border: none;
            background: #ffffff;
        }
        .modal__link {
            text-decoration: none;
            color: #3182ce;
            font-weight: 600;
        }
    """

    report_style_block = build_style_block(base_css, modal_css)

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

        if (!figure || !modal) {
            return;
        }

        const formatPercent = (value) => {
            if (value === null || value === undefined || value === "") {
                return "";
            }
            const text = String(value);
            return text.includes("%") ? text : `${text}%`;
        };

        const buildTitle = (point) => {
            const custom = Array.isArray(point.customdata) ? point.customdata : [];
            const label = custom[1];
            const allele = custom[2];
            const type = custom[3];
            const percent = formatPercent(custom[4]);
            const parts = [point.x, label, allele, type, percent].filter(Boolean);
            return parts.join(" ");
        };

        const openModal = (path, title) => {
            modal.classList.add("is-open");
            modal.setAttribute("aria-hidden", "false");
            modalFrame.src = path;
            modalTitle.textContent = title || "Allele report";
            modalLink.href = path;
        };

        const closeModal = () => {
            modal.classList.remove("is-open");
            modal.setAttribute("aria-hidden", "true");
            modalFrame.src = "";
            modalTitle.textContent = "Allele report";
            modalLink.removeAttribute("href");
        };

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

    inject_plot_assets(
        report_html_path,
        report_style_block,
        [report_hint_block, controls_block, modal_block],
        [script_block, report_script_block],
    )


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
    fileio.write_xlsx(results_summary, Path(report_directory, "read_summary.xlsx"))

    # Write to plot as HTML and PDF
    output_plot(results_summary, report_directory)
