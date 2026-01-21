(function() {
        const figure = document.getElementById("__PLOT_DIV_ID__");
        if (!figure) {
            return;
        }

        const fontInput = document.getElementById("plot-font-size");
        const applyButton = document.getElementById("plot-style-apply");
        const resetButton = document.getElementById("plot-style-reset");
        const colorPickerList = document.getElementById("color-picker-list");
        const exportControls = document.getElementById("export-controls");

        const clone = (obj) => JSON.parse(JSON.stringify(obj ?? (Array.isArray(obj) ? [] : {})));
        const initialData = clone(figure.data);
        const initialLayout = clone(figure.layout);
        const initialFontValue = fontInput.value;

        const traceGroups = __TRACE_GROUPS__;
        const exportButtons = __EXPORT_BUTTONS__;

        const ensureHex = (color) => {
            if (typeof color !== "string") {
                return "#1f77b4";
            }
            if (color.startsWith("#")) {
                if (color.length === 4) {
                    return "#" + color[1] + color[1] + color[2] + color[2] + color[3] + color[3];
                }
                return color;
            }
            const rgbMatch = color.match(/rgba?\\(([^)]+)\\)/i);
            if (rgbMatch) {
                const parts = rgbMatch[1].split(",").map((item) => parseFloat(item.trim()));
                const toHex = (value) => {
                    const v = Math.max(0, Math.min(255, Math.round(value)));
                    return v.toString(16).padStart(2, "0");
                };
                return "#" + toHex(parts[0]) + toHex(parts[1]) + toHex(parts[2]);
            }
            return "#1f77b4";
        };

        const getTraceColor = (trace) => {
            if (!trace || !trace.marker) {
                return "#1f77b4";
            }
            const color = trace.marker.color;
            if (Array.isArray(color) && color.length > 0) {
                return ensureHex(color[0]);
            }
            if (typeof color === "string") {
                return ensureHex(color);
            }
            return "#1f77b4";
        };

        const buildColorPickers = () => {
            colorPickerList.innerHTML = "";
            Object.entries(traceGroups).forEach(([name, indices]) => {
                const trace = figure.data[indices[0]];

                const wrapper = document.createElement("label");
                wrapper.className = "color-picker-item";

                const title = document.createElement("span");
                title.textContent = name;

                const picker = document.createElement("input");
                picker.type = "color";
                picker.value = getTraceColor(trace);

                picker.addEventListener("input", () => {
                    Plotly.restyle(figure, {"marker.color": picker.value}, indices);
                });
                title.addEventListener("click", () => picker.click());

                wrapper.appendChild(title);
                wrapper.appendChild(picker);
                colorPickerList.appendChild(wrapper);
            });
        };

        const applyFontSize = () => {
            const fontSize = parseInt(fontInput.value, 10);
            if (Number.isNaN(fontSize) || fontSize <= 0) {
                return;
            }
            Plotly.relayout(figure, {
                "font.size": fontSize,
                "xaxis.title.font.size": fontSize,
                "xaxis.tickfont.size": fontSize,
                "yaxis.title.font.size": fontSize,
                "yaxis.tickfont.size": fontSize,
                "legend.font.size": fontSize
            });
            Plotly.restyle(figure, {"textfont.size": fontSize});
        };

        const resetStyle = () => {
            Plotly.react(figure, clone(initialData), clone(initialLayout));
            fontInput.value = initialFontValue;
            buildColorPickers();
        };

        const buildExportButtons = () => {
            exportControls.innerHTML = "";
            exportButtons.forEach((btn) => {
                const button = document.createElement("button");
                button.textContent = btn.label;
                button.dataset.format = btn.format;
                button.addEventListener("click", async () => {
                    try {
                        await Plotly.downloadImage(figure, {
                            format: btn.format,
                            filename: "read_plot",
                            scale: 600 / 96,
                        });
                    } catch (error) {
                        console.error("Failed to export image", error);
                        alert("Failed to export " + btn.format.toUpperCase() + " file. Please ensure Plotly is fully loaded.");
                    }
                });
                exportControls.appendChild(button);
            });
        };

        applyButton.addEventListener("click", applyFontSize);
        resetButton.addEventListener("click", resetStyle);
        buildColorPickers();
        buildExportButtons();
    })();
