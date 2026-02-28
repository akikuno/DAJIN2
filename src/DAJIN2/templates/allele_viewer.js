(function() {
        const params = new URLSearchParams(window.location.search);
        const htmlPath = params.get("html") || "";
        const bamPath = params.get("bam") || "";
        const baiPath = params.get("bai") || "";
        const vcfPath = params.get("vcf") || "";
        const sample = params.get("sample") || "";
        const title = params.get("title") || "Allele report";
        const stateKey = params.get("state_key") || "";

        const frame = document.getElementById("viewer-frame");
        const status = document.getElementById("viewer-status");
        let igvView = document.getElementById("viewer-igv-view");
        const divider = document.getElementById("viewer-divider");
        const igvContainer = document.getElementById("viewer-igv");
        const viewer = document.getElementById("viewer");
        const viewerTitle = document.getElementById("viewer-title");
        const toolbar = document.getElementById("viewer-toolbar");
        const alleleSelect = document.getElementById("viewer-allele-select");
        const selectorNote = document.getElementById("viewer-selector-note");

        let splitRatio = 0.6;
        const assetPrefix = __ASSET_PREFIX__;
        const igvHelpers = window.DAJIN2IgvHelpers;

        if (!igvHelpers) {
            if (status) {
                status.textContent = "IGV helper is not available.";
            }
            return;
        }

        const resolvePath = (path) => igvHelpers.resolvePath(path, assetPrefix);
        const safeJsonParse = (text, fallback = null) => {
            try {
                return JSON.parse(text);
            } catch (_error) {
                return fallback;
            }
        };

        const normalizePercent = (value) => {
            const text = String(value ?? "").trim();
            if (!text) {
                return "";
            }
            return text.endsWith("%") ? text : `${text}%`;
        };

        const buildIgvPathsFromHtml = (path, sampleHint = "") => {
            if (!path) {
                return null;
            }
            const parts = String(path).split("/");
            const htmlIndex = parts.indexOf("HTML");
            const baseIndex = htmlIndex >= 0 ? htmlIndex : 0;
            if (parts.length <= baseIndex + 2) {
                return null;
            }
            const sampleName = parts[baseIndex + 1] || sampleHint;
            const filenameRaw = parts[parts.length - 1] || "";
            const filename = filenameRaw.split("?")[0].split("#")[0].trim();
            const stem = filename.replace(/\.html?$/i, "");
            const prefix = `${sampleName}_`;
            const header = stem.startsWith(prefix) ? stem.slice(prefix.length) : stem;
            if (!sampleName || !header) {
                return null;
            }
            return {
                sample: sampleName,
                bam: [".igvjs", sampleName, `${header}.bam`].join("/"),
                bai: [".igvjs", sampleName, `${header}.bam.bai`].join("/"),
                vcf: ["VCF", sampleName, `${sampleName}_${header}.vcf`].join("/"),
            };
        };

        const normalizeMember = (member, fallbackSample, fallbackTitle) => {
            const html = member?.path ? String(member.path) : member?.html ? String(member.html) : "";
            const igvPaths = buildIgvPathsFromHtml(html, fallbackSample);
            return {
                html,
                bam: member?.bam ? String(member.bam) : igvPaths?.bam || "",
                bai: member?.bai ? String(member.bai) : igvPaths?.bai || "",
                vcf: member?.vcf ? String(member.vcf) : igvPaths?.vcf || "",
                sample: member?.sample ? String(member.sample) : igvPaths?.sample || fallbackSample,
                title: member?.title ? String(member.title) : fallbackTitle,
                percent: normalizePercent(member?.percent),
            };
        };

        const getEntriesFromState = () => {
            if (!stateKey) {
                return [];
            }
            try {
                const text = window.localStorage.getItem(stateKey);
                if (!text) {
                    return [];
                }
                const state = safeJsonParse(text, {});
                if (Array.isArray(state?.members)) {
                    const fallbackTitle = state?.title ? String(state.title) : title;
                    return state.members
                        .map((member) => normalizeMember(member, sample, fallbackTitle))
                        .filter((entry) => entry.html);
                }
                return [];
            } catch (_error) {
                return [];
            }
        };

        const getEntriesFromQuery = () => {
            const payload = params.get("members") || "";
            if (!payload) {
                return [];
            }
            const parsed = safeJsonParse(payload, []);
            if (!Array.isArray(parsed)) {
                return [];
            }
            return parsed.map((member) => normalizeMember(member, sample, title)).filter((entry) => entry.html);
        };

        const buildFallbackEntry = () => {
            if (!htmlPath) {
                return null;
            }
            return {
                html: htmlPath,
                bam: bamPath,
                bai: baiPath,
                vcf: vcfPath,
                sample,
                title,
                percent: "",
            };
        };

        const dedupeEntries = (entries) => {
            const seen = new Set();
            return entries.filter((entry) => {
                const key = [entry.html, entry.bam, entry.bai, entry.vcf, entry.title].join("|");
                if (seen.has(key)) {
                    return false;
                }
                seen.add(key);
                return true;
            });
        };

        const stateEntries = getEntriesFromState();
        const queryEntries = getEntriesFromQuery();
        const fallbackEntry = buildFallbackEntry();
        const entries = dedupeEntries([
            ...stateEntries,
            ...queryEntries,
            ...(fallbackEntry ? [fallbackEntry] : []),
        ]);

        const applySplit = (topHeight) => {
            if (!viewer || !frame || !igvContainer || !divider) {
                return;
            }
            const toolbarHeight = toolbar && !toolbar.hidden ? toolbar.offsetHeight : 0;
            const total = viewer.clientHeight - divider.offsetHeight - toolbarHeight;
            if (total <= 0) {
                return;
            }
            const minTop = 180;
            const minBottom = 180;
            const maxTop = Math.max(minTop, total - minBottom);
            const preferredTop = total < minTop + minBottom ? total * splitRatio : topHeight;
            const upperBound = Math.max(0, Math.min(maxTop, total));
            const lowerBound = Math.max(0, Math.min(minTop, upperBound));
            const clampedTop = Math.max(lowerBound, Math.min(upperBound, preferredTop));
            frame.style.flex = "0 0 auto";
            igvContainer.style.flex = "0 0 auto";
            frame.style.height = `${clampedTop}px`;
            igvContainer.style.height = `${Math.max(0, total - clampedTop)}px`;
            splitRatio = total > 0 ? clampedTop / total : splitRatio;
        };

        if (divider && viewer) {
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
                if (event && divider.hasPointerCapture?.(event.pointerId)) {
                    divider.releasePointerCapture(event.pointerId);
                }
            };
            divider.addEventListener("pointerdown", (event) => {
                isDragging = true;
                startY = event.clientY;
                startTop = frame ? frame.getBoundingClientRect().height : 0;
                document.body.style.cursor = "row-resize";
                event.preventDefault();
                if (divider.setPointerCapture) {
                    divider.setPointerCapture(event.pointerId);
                }
            });
            divider.addEventListener("pointermove", onMove);
            divider.addEventListener("pointerup", stopDragging);
            divider.addEventListener("pointercancel", stopDragging);
            window.addEventListener("resize", () => {
                if (viewer.clientHeight > 0) {
                    applySplit(viewer.clientHeight * splitRatio);
                }
            });
        }

        const resetIgvView = () => {
            if (!igvContainer) {
                return;
            }
            if (igvView && igvView.parentNode) {
                igvView.parentNode.removeChild(igvView);
            }
            igvView = document.createElement("div");
            igvView.id = "viewer-igv-view";
            igvView.className = "viewer__igv-view";
            igvContainer.appendChild(igvView);
        };

        let igvBrowser = null;
        let igvRequestId = 0;

        const buildIgvBaseOptions = (sampleName) => {
            const igvGenomeInfo = __GENOME_INFO__;
            const baseOptions = igvGenomeInfo
                ? { genome: igvGenomeInfo.genome, locus: igvGenomeInfo.locus }
                : {
                      reference: {
                          fastaURL: sampleName ? resolvePath(`FASTA/${sampleName}/control.fasta`) : "",
                          indexURL: sampleName ? resolvePath(`FASTA/${sampleName}/control.fasta.fai`) : "",
                      },
                  };
            return { ...baseOptions, tracks: [] };
        };

        const loadIgv = async (entry) => {
            if (!status || !igvView) {
                return;
            }
            const hasBam = Boolean(entry.bam && entry.bai);
            const hasVcf = Boolean(entry.vcf);
            if (!hasBam && !hasVcf) {
                status.textContent = "No IGV tracks available.";
                igvView.innerHTML = "";
                return;
            }
            if (!window.igv) {
                status.textContent = "IGV library is not available.";
                igvView.innerHTML = "";
                return;
            }

            status.textContent = "Loading IGV...";
            const options = buildIgvBaseOptions(entry.sample);
            const currentRequestId = ++igvRequestId;
            try {
                const alignmentTrack = hasBam
                    ? igvHelpers.buildAlignmentTrack(resolvePath(entry.bam), resolvePath(entry.bai))
                    : null;
                const variantTrack = hasVcf ? igvHelpers.buildVcfTrack(resolvePath(entry.vcf)) : null;

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

                if (browser && typeof browser.loadTrack === "function") {
                    await igvHelpers.loadTracksVcfThenAlignment(browser, variantTrack, alignmentTrack);
                }

                igvBrowser = browser;
                status.textContent = "";
            } catch (error) {
                console.log("[DAJIN2][IGV][viewer] Failed to load IGV", {
                    error,
                    entry,
                    options,
                });
                status.textContent = "Failed to load IGV track.";
            }
        };

        const renderEntry = (entry, index, totalCount) => {
            if (frame) {
                frame.src = entry.html ? resolvePath(entry.html) : "";
            }
            if (viewerTitle) {
                const currentTitle = entry.title || title;
                viewerTitle.textContent = totalCount > 1 ? `${currentTitle} (${index + 1}/${totalCount})` : currentTitle;
            }
            loadIgv(entry);
            requestAnimationFrame(() => {
                if (viewer && viewer.clientHeight > 0) {
                    applySplit(viewer.clientHeight * splitRatio);
                }
            });
        };

        if (!entries.length) {
            if (status) {
                status.textContent = "No detailed allele report available.";
            }
            if (viewerTitle) {
                viewerTitle.textContent = title;
            }
            return;
        }

        if (toolbar && alleleSelect && selectorNote) {
            if (entries.length > 1) {
                toolbar.hidden = false;
                selectorNote.textContent = `Merged group with ${entries.length} alleles.`;
                alleleSelect.innerHTML = "";
                entries.forEach((entry, index) => {
                    const option = document.createElement("option");
                    option.value = String(index);
                    option.textContent = entry.percent
                        ? `${entry.title || `Allele ${index + 1}`} (${entry.percent})`
                        : entry.title || `Allele ${index + 1}`;
                    alleleSelect.appendChild(option);
                });
                alleleSelect.addEventListener("change", () => {
                    const index = Number.parseInt(alleleSelect.value, 10);
                    const safeIndex = Number.isFinite(index) ? index : 0;
                    const selected = entries[safeIndex] || entries[0];
                    renderEntry(selected, safeIndex, entries.length);
                });
            } else {
                toolbar.hidden = true;
            }
        }

        renderEntry(entries[0], 0, entries.length);
    })();
