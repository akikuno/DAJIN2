
def replace_N_to_match(cssplit, sequence):
    for i, cs in enumerate(cssplit):
        if not cs.startswith("N"):
            break
        cssplit[i] = "=" + sequence[i]
    cssplit = cssplit[::-1]
    for i, cs in enumerate(cssplit):
        if not cs.startswith("N"):
            break
        cssplit[i] = "=" + sequence[::-1][i]
    cssplit = cssplit[::-1]
    return cssplit


def call_cssplit_diffloci(clust_sample, DIFFLOCI_ALLELES, MASKS_CONTROL, DICT_ALLELE):
    for cssplit in clust_sample:
        qname = cssplit["QNAME"]
        allele = cssplit["ALLELE"]
        sv = cssplit["SV"]
        cssplit = cssplit["CSSPLIT"].split(",")
        diffloci = [loci for loci, _ in DIFFLOCI_ALLELES[f"{allele}-{sv}"]]
        maskloci = MASKS_CONTROL[allele]
        sequence = DICT_ALLELE[allele]
        if maskloci:
            cssplit = ["N" if mask else cs for cs, mask in zip(cssplit, maskloci)]
        cssplit = replace_N_to_match(cssplit, sequence)
