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
        const assetPrefix = __ASSET_PREFIX__;

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

        const normalizePercent = (value) => {
            if (value === null || value === undefined || value === "") {
                return "";
            }
            const text = String(value).trim();
            if (!text) {
                return "";
            }
            return text.endsWith("%") ? text : `${text}%`;
        };

        const buildHeaderFromCustom = (custom) => {
            if (!Array.isArray(custom) || custom.length < 5) {
                return "";
            }
            const label = custom[1];
            const allele = custom[2];
            const type = custom[3];
            const percent = custom[4];
            const labelText = label ? String(label).toLowerCase() : "";
            const alleleText = allele ? String(allele) : "";
            const typeText = type ? String(type) : "";
            const percentText = normalizePercent(percent);
            const header = [labelText, alleleText, typeText, percentText].filter(Boolean).join("_");
            return header.replace(/\s+/g, "");
        };

        const buildIgvPaths = (path, custom = null) => {
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
            const filenameRaw = parts[parts.length - 1];
            const filename = filenameRaw.split("?")[0].split("#")[0].trim();
            const stem = filename.replace(/\\.html?$/i, "");
            const prefix = `${sample}_`;
            const headerFromFile = stem.startsWith(prefix) ? stem.slice(prefix.length) : stem;
            const headerFromCustom = buildHeaderFromCustom(custom);
            const header = headerFromCustom || headerFromFile;
            return {
                sample,
                bam: [".igvjs", sample, `${header}.bam`].join("/"),
                bai: [".igvjs", sample, `${header}.bam.bai`].join("/"),
                vcf: ["VCF", sample, `${sample}_${header}.vcf`].join("/"),
            };
        };

        const buildViewerUrl = (path, title, custom = null) => {
            const igvPaths = buildIgvPaths(path, custom);
            if (!igvPaths) {
                return resolvePath(path);
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
                          fastaURL: resolvePath(`FASTA/${sample}/control.fasta`),
                          indexURL: resolvePath(`FASTA/${sample}/control.fasta.fai`),
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
            return rawItems
                .filter((item) => item.value !== undefined && item.value !== null && item.value !== "")
                .map((item) => ({ name: item.name, value: formatVcfValue(item.value, item.name) }));
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

        const buildVcfTrack = (vcfUrl, trackName) => ({
            name: trackName ? `${trackName} variants` : "Variants",
            url: vcfUrl,
            type: "variant",
            format: "vcf",
            displayMode: "EXPANDED",
            showAllSites: true,
            popupData: buildVcfPopupData,
            color: getVariantColor,
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

        const updateIgv = async (path, title, custom = null) => {
            if (!igvView || !igvStatus) {
                return;
            }
            const igvPaths = buildIgvPaths(path, custom);
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
                    const encodedTrackUrl = addCacheBuster(resolvePath(track.url), currentRequestId);
                    const encodedIndexUrl = addCacheBuster(resolvePath(track.indexURL), currentRequestId);
                    track.url = encodedTrackUrl;
                    track.indexURL = encodedIndexUrl;
                    track.indexurl = encodedIndexUrl;
                }
                const variantTrack = hasVcf ? buildVcfTrack(igvPaths.vcf, title) : null;
                if (variantTrack) {
                    variantTrack.url = addCacheBuster(resolvePath(variantTrack.url), currentRequestId);
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

        const openModal = (path, title, custom = null) => {
            modal.classList.add("is-open");
            modal.setAttribute("aria-hidden", "false");
            const encodedPath = resolvePath(path);
            modalFrame.src = encodedPath;
            modalTitle.textContent = title || "Allele report";
            modalLink.href = buildViewerUrl(path, title, custom);
            updateIgv(path, title, custom);
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

        const handlePlotlyClick = (event) => {
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
            openModal(detailPath, buildTitle(point), custom);
        };

        const attachPlotlyHandler = () => {
            if (figure && typeof figure.on === "function") {
                figure.on("plotly_click", handlePlotlyClick);
                return true;
            }
            return false;
        };

        const waitForPlotly = (timeoutMs = 5000, intervalMs = 100) =>
            new Promise((resolve) => {
                const start = Date.now();
                const tick = () => {
                    if (attachPlotlyHandler()) {
                        resolve(true);
                        return;
                    }
                if (Date.now() - start >= timeoutMs) {
                    resolve(false);
                    return;
                }
                    setTimeout(tick, intervalMs);
                };
                tick();
            });

        waitForPlotly();
    })();
