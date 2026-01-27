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
        const igvHelpers = window.DAJIN2IgvHelpers;
        const resolvePath = (path) =>
            igvHelpers && typeof igvHelpers.resolvePath === "function"
                ? igvHelpers.resolvePath(path, assetPrefix)
                : path;

        if (!figure || !modal) {
            return;
        }

        const formatPercent = (value) => {
            if (value === null || value === undefined || value === "") {
                return "";
            }
            return `${value}%`;
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

        const buildAlignmentTrack = (bamUrl, baiUrl, trackName) => {
            if (igvHelpers && typeof igvHelpers.buildAlignmentTrack === "function") {
                return igvHelpers.buildAlignmentTrack(bamUrl, baiUrl, trackName);
            }
            return null;
        };

        const buildVcfTrack = (vcfUrl, trackName) => {
            if (igvHelpers && typeof igvHelpers.buildVcfTrack === "function") {
                return igvHelpers.buildVcfTrack(vcfUrl, trackName);
            }
            return null;
        };

        const addCacheBuster = (url, seed) => {
            if (igvHelpers && typeof igvHelpers.addCacheBuster === "function") {
                return igvHelpers.addCacheBuster(url, seed);
            }
            if (!url) {
                return url;
            }
            const joiner = url.includes("?") ? "&" : "?";
            return `${url}${joiner}v=${seed}`;
        };

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
                const alignmentTrack = hasBam ? buildAlignmentTrack(igvPaths.bam, igvPaths.bai, title) : null;
                if (alignmentTrack) {
                    const encodedTrackUrl = addCacheBuster(resolvePath(alignmentTrack.url), currentRequestId);
                    const encodedIndexUrl = addCacheBuster(resolvePath(alignmentTrack.indexURL), currentRequestId);
                    alignmentTrack.url = encodedTrackUrl;
                    alignmentTrack.indexURL = encodedIndexUrl;
                    alignmentTrack.indexurl = encodedIndexUrl;
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
                if (browser && typeof browser.loadTrack === "function") {
                    if (igvHelpers && typeof igvHelpers.loadTracksVcfThenAlignment === "function") {
                        await igvHelpers.loadTracksVcfThenAlignment(browser, variantTrack, alignmentTrack);
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
