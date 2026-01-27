(function() {
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
        const assetPrefix = __ASSET_PREFIX__;

        const encodePath = (path) => {
            if (!path) {
                return "";
            }
            return path
                .split("/")
                .map((segment) => encodeURIComponent(segment))
                .join("/");
        };
        const resolvePath = (path) => {
            if (!path) {
                return "";
            }
            if (/^[a-zA-Z]+:\/\//.test(path) || path.startsWith("/")) {
                return path;
            }
            const combined = assetPrefix ? [assetPrefix, path].join("/") : path;
            return encodePath(combined);
        };

        if (frame) {
            frame.src = resolvePath(htmlPath);
        }
        if (document.getElementById("viewer-title")) {
            document.getElementById("viewer-title").textContent = title;
        }

        const applySplit = (topHeight) => {
            if (!viewer || !frame || !igvContainer || !divider) {
                return;
            }
            const total = viewer.clientHeight - divider.offsetHeight;
            const minTop = 180;
            const minBottom = 180;
            const clampedTop = Math.max(minTop, Math.min(total - minBottom, topHeight));
            frame.style.flex = "0 0 auto";
            igvContainer.style.flex = "0 0 auto";
            frame.style.height = `${clampedTop}px`;
            igvContainer.style.height = `${total - clampedTop}px`;
            splitRatio = clampedTop / total;
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
            if (viewer.clientHeight > 0) {
                requestAnimationFrame(() => {
                    applySplit(viewer.clientHeight * splitRatio);
                });
            }
        }

        const hasBam = Boolean(bamPath && baiPath);
        const hasVcf = Boolean(vcfPath);
        if (!hasBam && !hasVcf) {
            if (status) {
                status.textContent = "No IGV tracks available.";
            }
            return;
        }
        if (!window.igv) {
            if (status) {
                status.textContent = "IGV library is not available.";
            }
            return;
        }

        const igvGenomeInfo = __GENOME_INFO__;
        const baseOptions = igvGenomeInfo
            ? { genome: igvGenomeInfo.genome, locus: igvGenomeInfo.locus }
            : {
                  reference: {
                      fastaURL: sample ? resolvePath(`FASTA/${sample}/control.fasta`) : "",
                      indexURL: sample ? resolvePath(`FASTA/${sample}/control.fasta.fai`) : "",
                  },
              };
        const options = { ...baseOptions, tracks: [] };
        const alignmentTrack = hasBam
            ? {
                  name: title,
                  url: resolvePath(bamPath),
                  indexURL: resolvePath(baiPath),
                  indexurl: resolvePath(baiPath),
                  type: "alignment",
                  format: "bam",
                  autoHeight: true,
                  viewAsPairs: true,
                  samplingDepth: 30,
                  showInsertionText: true,
                  showDeletionText: true,
              }
            : null;
        const formatVcfValue = (value, context) => {
            if (value === null || value === undefined) {
                return "";
            }
            const text = Array.isArray(value)
                ? value.join(",")
                : typeof value === "object"
                  ? JSON.stringify(value)
                  : String(value);
            return text;
        };

        const buildPopupItems = (rawItems) => {
            const items = rawItems
                .filter((item) => item.value !== undefined && item.value !== null && item.value !== "")
                .map((item) => ({
                    name: item.name,
                    value: formatVcfValue(item.value, item.name),
                    rawValue: item.value,
                }));
            return items.map(({ name, value }) => ({ name, value }));
        };

        const normalizeVariantType = (value) => {
            if (value === null || value === undefined) {
                return "";
            }
            return String(value).toUpperCase();
        };

        const resolveVariantType = (feature) => {
            if (!feature) {
                return "";
            }
            const info = feature.info || feature.INFO || {};
            const svtype = normalizeVariantType(info.SVTYPE || feature.SVTYPE || feature.svtype);
            const type = normalizeVariantType(info.TYPE || feature.TYPE || feature.type);
            if (svtype) {
                return svtype;
            }
            if (type) {
                return type;
            }
            const alt = normalizeVariantType(feature.alt || feature.ALT || "");
            if (alt.includes("<INV>")) {
                return "INV";
            }
            if (alt.includes("<DEL>")) {
                return "DEL";
            }
            if (alt.includes("<INS>")) {
                return "INS";
            }
            return "";
        };

        const getVariantColor = (feature) => {
            const colors = {
                SUB: "#2ca02c",
                INS: "#d62728",
                DEL: "#1f77b4",
                INV: "#9467bd",
            };
            return colors[resolveVariantType(feature)] || "#4a5568";
        };

        const extractInfoEntries = (info) => {
            if (!info) {
                return [];
            }
            if (info instanceof Map) {
                const entries = [];
                info.forEach((value, key) => {
                    entries.push([key, value]);
                });
                return entries;
            }
            if (typeof info === "object") {
                return Object.entries(info);
            }
            return [];
        };

        const extractAltValue = (feature) => {
            return (
                feature.alt ||
                feature.ALT ||
                feature.altBases ||
                feature.alts ||
                feature.alleles ||
                feature.altBasesString ||
                ""
            );
        };

        const buildVcfPopupData = (feature) => {
            if (!feature) {
                return [{ name: "Variant", value: "No data" }];
            }
            const chr = feature.chr || feature.chrom || feature.CHROM || "";
            const pos =
                feature.pos ?? feature.POS ?? (feature.start != null ? feature.start + 1 : feature.position ?? "");
            const ref = feature.ref || feature.REF || "";
            const alt = extractAltValue(feature);
            const id = feature.id || feature.ID || "";
            const qual = feature.qual || feature.QUAL || "";
            const filter = feature.filter || feature.FILTER || "";
            const info = feature.info || feature.INFO || {};
            const items = [];
            if (chr) {
                items.push({ name: "CHROM", value: chr });
            }
            if (pos !== "") {
                items.push({ name: "POS", value: pos });
            }
            if (ref) {
                items.push({ name: "REF", value: ref });
            }
            if (alt) {
                items.push({ name: "ALT", value: alt });
            }
            if (id && id !== ".") {
                items.push({ name: "ID", value: id });
            }
            if (qual !== "" && qual !== ".") {
                items.push({ name: "QUAL", value: qual });
            }
            if (filter && filter !== ".") {
                items.push({ name: "FILTER", value: filter });
            }
            extractInfoEntries(info)
                .sort(([a], [b]) => String(a).localeCompare(String(b)))
                .forEach(([key, value]) => {
                    items.push({ name: String(key), value });
                });
            if (items.length === 0) {
                return buildPopupItems([{ name: "Variant", value: "No details" }]);
            }
            return buildPopupItems(items);
        };

        const variantTrack = hasVcf
            ? {
                  name: `${title} variants`,
                  url: resolvePath(vcfPath),
                  type: "variant",
                  format: "vcf",
                  displayMode: "EXPANDED",
                  showAllSites: true,
                  popupData: buildVcfPopupData,
                  color: getVariantColor,
              }
            : null;

        igv.createBrowser(igvView, options)
            .then((browser) => {
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
                const enforceTrackOrder = async (target, trackConfig) => {
                    if (!target || !Array.isArray(target.tracks)) {
                        return;
                    }
                    const refseqIndex = target.tracks.findIndex((t) =>
                        String(t?.name ?? "").toLowerCase().includes("refseq curated")
                    );
                    const alignmentTrackLoaded = target.tracks.find(
                        (t) => t?.type === "alignment" && String(t?.name ?? "") === String(trackConfig?.name ?? "")
                    );
                    if (refseqIndex === -1 || !alignmentTrackLoaded) {
                        return;
                    }
                    const alignmentIndex = target.tracks.indexOf(alignmentTrackLoaded);
                    if (alignmentIndex > refseqIndex) {
                        return;
                    }
                    if (typeof target.removeTrack === "function") {
                        target.removeTrack(alignmentTrackLoaded);
                        if (typeof target.loadTrack === "function") {
                            await target.loadTrack(trackConfig);
                        }
                    }
                };
                const loadTrackSafe = async (target, trackConfig) => {
                    if (!trackConfig || typeof target?.loadTrack !== "function") {
                        return;
                    }
                    await target.loadTrack(trackConfig);
                };
                if (browser && typeof browser.loadTrack === "function") {
                    return waitForRefseq(browser)
                        .then(() => loadTrackSafe(browser, alignmentTrack))
                        .then(() => (alignmentTrack ? enforceTrackOrder(browser, alignmentTrack) : null))
                        .then(() => loadTrackSafe(browser, variantTrack))
                        .then(() => browser);
                }
                return browser;
            })
            .then((browser) => {
                if (status) {
                    status.textContent = "";
                }
            })
            .catch((error) => {
                console.error("Failed to load IGV", error);
                if (status) {
                    status.textContent = "Failed to load IGV track.";
                }
            });
    })();
