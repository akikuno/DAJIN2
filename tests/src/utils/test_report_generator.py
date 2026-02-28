from __future__ import annotations

import json
import subprocess
import textwrap
from pathlib import Path


def test_report_controls_deduplicate_same_allele_type_without_collapsing_stacks(tmp_path):
    template_path = Path(__file__).resolve().parents[3] / "src" / "DAJIN2" / "templates" / "report_controls.js"
    script = template_path.read_text(encoding="utf-8").replace("__PLOT_DIV_ID__", "read_plot_fig").replace(
        "__EXPORT_BUTTONS__",
        "[]",
    )

    initial_data = [
        {
            "name": "deletion with indels",
            "marker": {"color": "#636efa"},
            "x": ["sample", "sample", "sample", "sample"],
            "y": [33.884, 33.402, 31.749, 0.964],
            "customdata": [
                ["HTML/sample/sample_allele01_deletion_indels_33.884%.html", "Allele01", "deletion", "indels", 33.884],
                ["HTML/sample/sample_allele02_deletion_indels_33.402%.html", "Allele02", "deletion", "indels", 33.402],
                ["HTML/sample/sample_allele03_deletion_indels_31.749%.html", "Allele03", "deletion", "indels", 31.749],
                ["HTML/sample/sample_allele04_deletion_indels_0.964%.html", "Allele04", "deletion", "indels", 0.964],
            ],
        }
    ]

    harness = textwrap.dedent(
        f"""
        const scriptText = {json.dumps(script)};
        const initialData = {json.dumps(initial_data)};

        class Element {{
            constructor(tagName = "div", id = "") {{
                this.tagName = String(tagName).toUpperCase();
                this.id = id;
                this.children = [];
                this.listeners = {{}};
                this.style = {{}};
                this.dataset = {{}};
                this.textContent = "";
                this.value = "";
                this.type = "";
                this.checked = false;
                this.className = "";
                this.parentNode = null;
                this._innerHTML = "";
            }}

            appendChild(child) {{
                child.parentNode = this;
                this.children.push(child);
                return child;
            }}

            addEventListener(name, listener) {{
                this.listeners[name] = listener;
            }}

            click() {{
                if (this.listeners.click) {{
                    this.listeners.click({{ target: this }});
                }}
            }}

            querySelectorAll(selector) {{
                const matches = [];
                const walk = (node) => {{
                    if (
                        selector === 'input[type="checkbox"]:checked' &&
                        node.tagName === "INPUT" &&
                        node.type === "checkbox" &&
                        node.checked
                    ) {{
                        matches.push(node);
                    }}
                    node.children.forEach(walk);
                }};
                this.children.forEach(walk);
                return matches;
            }}

            set innerHTML(value) {{
                this._innerHTML = value;
                if (value === "") {{
                    this.children = [];
                }}
            }}

            get innerHTML() {{
                return this._innerHTML;
            }}
        }}

        const elements = new Map();
        const register = (id, element) => {{
            elements.set(id, element);
            return element;
        }};

        const figure = register("read_plot_fig", new Element("div", "read_plot_fig"));
        figure.data = initialData;
        figure.layout = {{ font: {{ size: 18 }} }};
        figure.on = () => {{}};

        register("plot-font-size", Object.assign(new Element("input", "plot-font-size"), {{ value: "18" }}));
        register("plot-style-apply", new Element("button", "plot-style-apply"));
        register("plot-style-reset", new Element("button", "plot-style-reset"));
        register("color-picker-list", new Element("div", "color-picker-list"));
        register("export-controls", new Element("div", "export-controls"));
        register("allele-group-list", new Element("div", "allele-group-list"));
        register("allele-rename-input", new Element("input", "allele-rename-input"));
        register("allele-rename-apply", new Element("button", "allele-rename-apply"));
        register("allele-editor-reset", new Element("button", "allele-editor-reset"));
        register("allele-editor-status", new Element("div", "allele-editor-status"));

        global.document = {{
            getElementById(id) {{
                return elements.get(id) || null;
            }},
            createElement(tagName) {{
                return new Element(tagName);
            }},
        }};

        global.window = {{}};
        global.Plotly = {{
            react(target, data, layout) {{
                global.__plotlyData = data;
                target.data = data;
                target.layout = layout;
            }},
            relayout() {{}},
            restyle() {{}},
            downloadImage() {{
                return Promise.resolve();
            }},
        }};

        eval(scriptText);

        const groupList = elements.get("allele-group-list");
        const colorList = elements.get("color-picker-list");
        const traces = global.__plotlyData || [];

        console.log(JSON.stringify({{
            groupCount: groupList.children.length,
            colorCount: colorList.children.length,
            traceCount: traces.length,
            legendCount: traces.filter((trace) => trace.showlegend).length,
            groupLabel: groupList.children[0]?.children[2]?.textContent || "",
            groupMeta: groupList.children[0]?.children[3]?.textContent || "",
        }}));
        """
    )

    harness_path = tmp_path / "report_controls_harness.js"
    harness_path.write_text(harness, encoding="utf-8")

    result = subprocess.run(["node", str(harness_path)], capture_output=True, text=True, check=False)

    assert result.returncode == 0, result.stderr
    output = json.loads(result.stdout)
    assert output["groupCount"] == 1
    assert output["colorCount"] == 1
    assert output["traceCount"] == 4
    assert output["legendCount"] == 1
    assert output["groupLabel"] == "deletion with indels"
    assert output["groupMeta"] == "4 alleles"
