(function () {
    const figure = document.getElementById("__PLOT_DIV_ID__");
    if (!figure) {
        return;
    }

    const fontInput = document.getElementById("plot-font-size");
    const applyButton = document.getElementById("plot-style-apply");
    const resetButton = document.getElementById("plot-style-reset");
    const colorPickerList = document.getElementById("color-picker-list");
    const exportControls = document.getElementById("export-controls");

    const alleleGroupList = document.getElementById("allele-group-list");
    const alleleRenameInput = document.getElementById("allele-rename-input");
    const alleleRenameApply = document.getElementById("allele-rename-apply");
    const alleleEditorReset = document.getElementById("allele-editor-reset");
    const alleleEditorStatus = document.getElementById("allele-editor-status");

    const clone = (obj) => JSON.parse(JSON.stringify(obj ?? (Array.isArray(obj) ? [] : {})));
    const initialData = clone(figure.data);
    const initialLayout = clone(figure.layout);
    const initialFontValue = fontInput ? fontInput.value : "18";

    const exportButtons = __EXPORT_BUTTONS__;

    const toNumber = (value) => {
        const number = Number(value);
        return Number.isFinite(number) ? number : 0;
    };

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
            const toHex = (channel) => {
                const valueClamped = Math.max(0, Math.min(255, Math.round(channel)));
                return valueClamped.toString(16).padStart(2, "0");
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

    const formatAlleleType = (allele, type) => {
        const alleleTextRaw = String(allele ?? "").trim();
        const typeTextRaw = String(type ?? "").trim();
        if (alleleTextRaw.includes("|")) {
            const parts = alleleTextRaw.split("|").map((part) => part.trim()).filter(Boolean);
            if (parts.length >= 2) {
                return parts[1].toLowerCase();
            }
            if (parts.length === 1) {
                return parts[0].toLowerCase();
            }
            return "";
        }
        const alleleText = alleleTextRaw;
        const typeText = typeTextRaw;
        const alleleLower = String(alleleText ?? "").trim().toLowerCase();
        const typeLower = String(typeText ?? "").trim().toLowerCase();
        if (!typeLower || typeLower === "intact") {
            return alleleLower;
        }
        if (typeLower === "indels") {
            return `${alleleLower} with indels`.trim();
        }
        return `${alleleLower} ${typeLower}`.trim();
    };

    const formatPercentLabel = (value) => {
        if (!Number.isFinite(value)) {
            return "0";
        }
        const rounded = Math.round(value * 100) / 100;
        return rounded.toFixed(2).replace(/\\.00$/, "").replace(/(\\.\\d)0$/, "$1");
    };

    const buildPointTitle = (sample, label, allele, type, percent) => {
        const alleleType = formatAlleleType(allele, type);
        const percentText = Number.isFinite(percent) ? `(${formatPercentLabel(percent)}%)` : "";
        const labelText = String(label ?? "").toLowerCase().trim();
        return [sample, alleleType || labelText, percentText].filter(Boolean).join(" ").replace(/\\s+/g, " ").trim();
    };

    const buildSourceModel = (plotData) => {
        const model = {};
        const displayColorMap = {};
        const buildSourceKey = (label, allele, type, traceName, fallbackIndex) => {
            const labelKey = String(label ?? "").trim().toLowerCase();
            const alleleKey = String(allele ?? "").trim().toLowerCase();
            const typeKey = String(type ?? "").trim().toLowerCase();
            const alleleDisplayKey = formatAlleleType(allele, type).toLowerCase();
            if (labelKey || alleleKey || typeKey) {
                return [labelKey || alleleKey, alleleDisplayKey, alleleKey, typeKey].join("|");
            }
            return `${String(traceName ?? "").trim().toLowerCase()}|${fallbackIndex}`;
        };

        const buildSourceName = (allele, type, traceName) => {
            const alleleType = formatAlleleType(allele, type);
            return alleleType || String(traceName ?? "").trim().toLowerCase();
        };

        (plotData || []).forEach((trace, traceIndex) => {
            const traceNameRaw = trace && trace.name ? String(trace.name) : `Series ${traceIndex + 1}`;
            const traceName = traceNameRaw.trim() || `Series ${traceIndex + 1}`;

            const xs = Array.isArray(trace?.x) ? trace.x : [];
            const ys = Array.isArray(trace?.y) ? trace.y : [];
            const customList = Array.isArray(trace?.customdata) ? trace.customdata : [];
            const length = Math.max(xs.length, ys.length);

            for (let index = 0; index < length; index += 1) {
                const custom = Array.isArray(customList[index]) ? customList[index] : [];
                const sample = String(xs[index] ?? "").trim();
                const value = toNumber(ys[index]);
                const label = custom[1] ? String(custom[1]).trim() : "";
                const allele = custom[2] ? String(custom[2]).trim() : traceName;
                const type = custom[3] ? String(custom[3]).trim() : "";
                const percentRaw = custom[4] !== undefined && custom[4] !== null && custom[4] !== ""
                    ? String(custom[4]).replace(/%$/, "")
                    : value;
                const percent = toNumber(percentRaw);
                const path = custom[0] ? String(custom[0]).trim() : "";
                const sourceKey = buildSourceKey(label, allele, type, traceName, index);
                if (!model[sourceKey]) {
                    const sourceName = buildSourceName(allele, type, traceName);
                    const colorKey = sourceName.toLowerCase();
                    if (!displayColorMap[colorKey]) {
                        displayColorMap[colorKey] = getTraceColor(trace);
                    }
                    model[sourceKey] = {
                        name: sourceName,
                        color: displayColorMap[colorKey],
                        points: [],
                    };
                }
                const title = buildPointTitle(sample, label, allele, type, percent);
                model[sourceKey].points.push({
                    sample,
                    value,
                    percent,
                    path,
                    label,
                    allele,
                    type,
                    title,
                    traceName: model[sourceKey].name,
                });
            }
        });
        return model;
    };

    const buildSampleOrder = (plotData, sourceModel) => {
        const order = [];
        const seen = new Set();
        (plotData || []).forEach((trace) => {
            const xs = Array.isArray(trace?.x) ? trace.x : [];
            xs.forEach((sample) => {
                const name = String(sample ?? "").trim();
                if (!name || seen.has(name)) {
                    return;
                }
                seen.add(name);
                order.push(name);
            });
        });
        if (order.length > 0) {
            return order;
        }
        Object.values(sourceModel).forEach((source) => {
            source.points.forEach((point) => {
                if (!point.sample || seen.has(point.sample)) {
                    return;
                }
                seen.add(point.sample);
                order.push(point.sample);
            });
        });
        return order;
    };

    const sourceModel = buildSourceModel(initialData);
    const sampleOrder = buildSampleOrder(initialData, sourceModel);
    const sourceNames = Object.keys(sourceModel);

    const normalizeGroupNameKey = (name) => String(name ?? "").trim().toLowerCase();
    const createGroupId = (index) => `group_${index + 1}`;

    const mergeGroupsByName = (inputGroups = []) => {
        const mergedGroups = [];
        const groupIndexByName = new Map();

        inputGroups.forEach((group, index) => {
            const members = Array.isArray(group?.members)
                ? group.members.map((member) => String(member ?? "").trim()).filter(Boolean)
                : [];
            if (!members.length) {
                return;
            }

            const fallbackName = sourceModel[members[0]]?.name || members[0];
            const name = String(group?.name ?? fallbackName).trim() || fallbackName;
            const color = String(group?.color ?? sourceModel[members[0]]?.color ?? "#1f77b4");
            const key = normalizeGroupNameKey(name) || createGroupId(index);

            if (!groupIndexByName.has(key)) {
                mergedGroups.push({
                    id: String(group?.id ?? createGroupId(mergedGroups.length)).trim() || createGroupId(mergedGroups.length),
                    name,
                    members: [...members],
                    color,
                });
                groupIndexByName.set(key, mergedGroups.length - 1);
                return;
            }

            const target = mergedGroups[groupIndexByName.get(key)];
            const memberSet = new Set(target.members);
            members.forEach((member) => {
                if (!memberSet.has(member)) {
                    memberSet.add(member);
                    target.members.push(member);
                }
            });
        });

        return mergedGroups;
    };

    const buildInitialGroups = () =>
        mergeGroupsByName(
            sourceNames.map((sourceKey, index) => ({
                id: createGroupId(index),
                name: sourceModel[sourceKey]?.name || sourceKey,
                members: [sourceKey],
                color: sourceModel[sourceKey]?.color || "#1f77b4",
            })),
        );

    const initialGroups = buildInitialGroups();
    let groups = clone(initialGroups);

    const syncGroupColorsByName = () => {
        groups = mergeGroupsByName(groups);
        const colorMap = {};
        groups.forEach((group) => {
            const key = normalizeGroupNameKey(group.name);
            if (!key) {
                return;
            }
            if (!colorMap[key]) {
                colorMap[key] = group.color;
            }
            group.color = colorMap[key];
        });
    };

    const getSelectedGroupIds = () => {
        if (!alleleGroupList) {
            return [];
        }
        return Array.from(alleleGroupList.querySelectorAll('input[type="checkbox"]:checked')).map((input) => String(input.value));
    };

    const setEditorStatus = (message, isError = false) => {
        if (!alleleEditorStatus) {
            return;
        }
        alleleEditorStatus.textContent = message || "";
        alleleEditorStatus.dataset.kind = isError ? "error" : "info";
    };

    const buildGroupList = (selectedIds = []) => {
        if (!alleleGroupList) {
            return;
        }
        const selectedSet = new Set(selectedIds);
        alleleGroupList.innerHTML = "";
        groups.forEach((group) => {
            const item = document.createElement("label");
            item.className = "allele-editor__item";

            const checkbox = document.createElement("input");
            checkbox.type = "checkbox";
            checkbox.value = group.id;
            checkbox.checked = selectedSet.has(group.id);

            const swatch = document.createElement("span");
            swatch.className = "allele-editor__swatch";
            swatch.style.backgroundColor = group.color;

            const name = document.createElement("span");
            name.className = "allele-editor__name";
            name.textContent = group.name;

            const meta = document.createElement("span");
            meta.className = "allele-editor__meta";
            const memberCount = Array.isArray(group.members) ? group.members.length : 0;
            meta.textContent = memberCount > 1 ? `${memberCount} alleles` : "";

            item.appendChild(checkbox);
            item.appendChild(swatch);
            item.appendChild(name);
            item.appendChild(meta);
            alleleGroupList.appendChild(item);
        });
    };

    const buildMembersPayload = (points) =>
        points.map((point) => ({
            path: point.path,
            label: point.label,
            allele: point.allele,
            type: point.type,
            percent: formatPercentLabel(point.percent),
            title: point.title,
            sample: point.sample,
            trace: point.traceName,
        }));

    const buildTooltipSummary = (points) => {
        if (!points.length) {
            return "No detailed allele report available.";
        }
        if (points.length === 1) {
            return `Detail: ${points[0].title}`;
        }
        const names = points.map((point) => formatAlleleType(point.allele, point.type) || point.traceName);
        const uniqueNames = Array.from(new Set(names));
        const preview = uniqueNames.slice(0, 3).join(", ");
        const extra = uniqueNames.length > 3 ? ` (+${uniqueNames.length - 3} more)` : "";
        return `Contains ${points.length} allele entries: ${preview}${extra}`;
    };

    const buildDisplayTraces = () => {
        const orderedSamples = sampleOrder.length > 0 ? sampleOrder : [];
        const traces = [];

        groups.forEach((group) => {
            group.members.forEach((memberName, memberIndex) => {
                const source = sourceModel[memberName];
                if (!source) {
                    return;
                }

                const buckets = new Map(orderedSamples.map((sample) => [sample, { total: 0, points: [] }]));
                source.points.forEach((point) => {
                    if (!point.sample) {
                        return;
                    }
                    if (!buckets.has(point.sample)) {
                        buckets.set(point.sample, { total: 0, points: [] });
                    }
                    const bucket = buckets.get(point.sample);
                    bucket.total += point.value;
                    if (point.path || point.value > 0) {
                        bucket.points.push(point);
                    }
                });

                const x = [];
                const y = [];
                const text = [];
                const customdata = [];

                Array.from(buckets.keys()).forEach((sample) => {
                    const bucket = buckets.get(sample) || { total: 0, points: [] };
                    const sortedPoints = bucket.points.slice().sort((a, b) => b.value - a.value);
                    const primary = sortedPoints.find((point) => point.path) || sortedPoints[0] || null;
                    const membersPayload = buildMembersPayload(sortedPoints);
                    const total = Math.round(bucket.total * 1000000) / 1000000;

                    x.push(sample);
                    y.push(total);
                    text.push(total > 0 ? formatPercentLabel(total) : "");
                    customdata.push([
                        primary?.path || "",
                        primary?.label || "",
                        primary?.allele || group.name,
                        primary?.type || "",
                        total,
                        JSON.stringify(membersPayload),
                        group.name,
                        buildTooltipSummary(sortedPoints),
                    ]);
                });

                traces.push({
                    type: "bar",
                    name: group.name,
                    legendgroup: group.id,
                    showlegend: memberIndex === 0,
                    marker: { color: group.color },
                    x,
                    y,
                    text,
                    textposition: "inside",
                    cliponaxis: false,
                    customdata,
                    hovertemplate: "Sample: %{x}<br>Allele: %{customdata[6]}<br>Percent: %{y:.2f}%<br>%{customdata[7]}<extra></extra>",
                });
            });
        });

        return traces;
    };

    const renderGroupedPlot = () => {
        syncGroupColorsByName();
        const layout = Object.keys(figure.layout || {}).length > 0 ? clone(figure.layout) : clone(initialLayout);
        Plotly.react(figure, buildDisplayTraces(), layout);
        buildColorPickers();
    };

    const buildColorPickers = () => {
        if (!colorPickerList) {
            return;
        }
        colorPickerList.innerHTML = "";
        groups.forEach((group) => {
            const wrapper = document.createElement("label");
            wrapper.className = "color-picker-item";

            const title = document.createElement("span");
            title.textContent = group.name;

            const picker = document.createElement("input");
            picker.type = "color";
            picker.value = ensureHex(group.color);

            picker.addEventListener("input", () => {
                const selectedIds = getSelectedGroupIds();
                const targetNameKey = normalizeGroupNameKey(group.name);
                groups.forEach((currentGroup) => {
                    if (normalizeGroupNameKey(currentGroup.name) === targetNameKey) {
                        currentGroup.color = picker.value;
                    }
                });
                renderGroupedPlot();
                buildGroupList(selectedIds);
            });
            title.addEventListener("click", () => picker.click());

            wrapper.appendChild(title);
            wrapper.appendChild(picker);
            colorPickerList.appendChild(wrapper);
        });
    };

    const applyFontSize = () => {
        if (!fontInput) {
            return;
        }
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
            "legend.font.size": fontSize,
        });
        Plotly.restyle(figure, { "textfont.size": fontSize });
    };

    const clearEditorInputs = () => {
        if (alleleRenameInput) {
            alleleRenameInput.value = "";
        }
    };

    const resetAlleleEdits = (message = "Allele edits were reset.") => {
        groups = clone(initialGroups);
        clearEditorInputs();
        buildGroupList();
        renderGroupedPlot();
        setEditorStatus(message);
    };

    const resetStyle = () => {
        Plotly.react(figure, clone(initialData), clone(initialLayout));
        if (fontInput) {
            fontInput.value = initialFontValue;
        }
        groups = clone(initialGroups);
        clearEditorInputs();
        buildGroupList();
        renderGroupedPlot();
        setEditorStatus("Styles and allele edits were reset.");
    };

    const applyRename = () => {
        const selectedIds = getSelectedGroupIds();
        if (selectedIds.length < 1) {
            setEditorStatus("Select one or more alleles to rename.", true);
            return;
        }
        const nextName = String(alleleRenameInput?.value ?? "").trim();
        if (!nextName) {
            setEditorStatus("Enter a new allele name.", true);
            return;
        }
        const targets = groups.filter((group) => selectedIds.includes(group.id));
        if (!targets.length) {
            setEditorStatus("Selected allele was not found.", true);
            return;
        }
        targets.forEach((group) => {
            group.name = nextName;
        });
        renderGroupedPlot();
        buildGroupList();
        setEditorStatus(`Renamed ${targets.length} alleles to \"${nextName}\".`);
    };

    const buildExportButtons = () => {
        if (!exportControls) {
            return;
        }
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

    applyButton?.addEventListener("click", applyFontSize);
    resetButton?.addEventListener("click", resetStyle);
    alleleRenameApply?.addEventListener("click", applyRename);
    alleleEditorReset?.addEventListener("click", () => resetAlleleEdits());

    buildGroupList();
    buildExportButtons();
    renderGroupedPlot();
})();
