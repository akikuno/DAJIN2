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
        const igvHelpers = window.DAJIN2IgvHelpers;
        if (!igvHelpers) {
            if (status) {
                status.textContent = "IGV helper is not available.";
            }
            return;
        }
        const resolvePath = (path) => igvHelpers.resolvePath(path, assetPrefix);

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
            ? igvHelpers.buildAlignmentTrack(resolvePath(bamPath), resolvePath(baiPath), title)
            : null;
        const variantTrack = hasVcf ? igvHelpers.buildVcfTrack(resolvePath(vcfPath), title) : null;

        igv.createBrowser(igvView, options)
            .then((browser) => {
                if (browser && typeof browser.loadTrack === "function") {
                    return igvHelpers.loadTracksVcfThenAlignment(browser, variantTrack, alignmentTrack).then(() => browser);
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
