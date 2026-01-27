(function() {
        const global = typeof window !== "undefined" ? window : globalThis;
        if (global.DAJIN2IgvHelpers) {
            return;
        }

        const encodePath = (path) => {
            if (!path) {
                return "";
            }
            return path
                .split("/")
                .map((segment) => encodeURIComponent(segment))
                .join("/");
        };

        const resolvePath = (path, assetPrefix) => {
            if (!path) {
                return "";
            }
            if (/^[a-zA-Z]+:\/\//.test(path) || path.startsWith("/")) {
                return path;
            }
            const combined = assetPrefix ? [assetPrefix, path].join("/") : path;
            return encodePath(combined);
        };

        const addCacheBuster = (url, seed) => {
            if (!url) {
                return url;
            }
            const joiner = url.includes("?") ? "&" : "?";
            return `${url}${joiner}v=${seed}`;
        };

        const buildAlignmentTrack = (bamUrl, baiUrl, trackName) => ({
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

        const formatVcfValue = (value) => {
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
                .map((item) => ({ name: item.name, value: formatVcfValue(item.value) }));
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
            popupData: (feature) => buildVcfPopupData(feature),
            color: getVariantColor,
        });

        const waitForRefseq = (target, timeoutMs = 2000) =>
            new Promise((resolve) => {
                const start = Date.now();
                const check = () => {
                    const hasRefseq = Array.isArray(target?.tracks)
                        ? target.tracks.some((t) => String(t?.name ?? "").toLowerCase().includes("refseq curated"))
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

        const enforceAlignmentAfterRefseq = async (browser, trackConfig) => {
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

        const loadTracksVcfThenAlignment = async (browser, variantTrack, alignmentTrack) => {
            if (!browser || typeof browser.loadTrack !== "function") {
                return;
            }
            await waitForRefseq(browser, 2000);
            if (variantTrack) {
                await browser.loadTrack(variantTrack);
            }
            if (alignmentTrack) {
                await browser.loadTrack(alignmentTrack);
                await enforceAlignmentAfterRefseq(browser, alignmentTrack);
            }
        };

        global.DAJIN2IgvHelpers = {
            encodePath,
            resolvePath,
            addCacheBuster,
            buildAlignmentTrack,
            buildVcfTrack,
            waitForRefseq,
            enforceAlignmentAfterRefseq,
            loadTracksVcfThenAlignment,
        };
})();
