#!--------------------------------------------------
print(f"marge_alleles.py: WOURKING IN PROGRESS...")
#!--------------------------------------------------


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


def replace_maskloci_to_match(cssplit, sequence, maskloci):
    for i, mask in enumerate(maskloci):
        if mask:
            cssplit[i] = "=" + sequence[i]
    return cssplit


def call_cssplit_diffloci(clust_sample, DIFFLOCI_ALLELES, MASKS_CONTROL, DICT_ALLELE):
    cssplit_diffloci = []
    for cssplit in clust_sample:
        qname = cssplit["QNAME"]
        label = cssplit["LABEL"]
        allele = cssplit["ALLELE"]
        sv = cssplit["SV"]
        cssplit = cssplit["CSSPLIT"].split(",")
        maskloci = MASKS_CONTROL[allele]
        diffloci = {loci for loci, _ in DIFFLOCI_ALLELES[f"{allele}-{sv}"]}
        sequence = DICT_ALLELE[allele]
        cssplit = replace_N_to_match(cssplit, sequence)
        if maskloci:
            cssplit = replace_maskloci_to_match(cssplit, sequence, maskloci)
        if diffloci:
            cssplit = [cs for i, cs in enumerate(cssplit) if i in diffloci]
        cssplit = "".join(cssplit)
        cssplit_diffloci.append({"QNAME": qname, "LABEL": label, "CSSPLIT_DIFFLOCI": cssplit})
    return cssplit_diffloci


# cssample = call_cssplit_diffloci(clust_sample, DIFFLOCI_ALLELES, MASKS_CONTROL, DICT_ALLELE)
# cscontrol = call_cssplit_diffloci(clust_control, DIFFLOCI_ALLELES, MASKS_CONTROL, DICT_ALLELE)

# cssample_set = {cs["CSSPLIT_DIFFLOCI"] for cs in cssample}
# cscontrol_set = {cs["CSSPLIT_DIFFLOCI"] for cs in cscontrol}

# len(cssample_set)

# len(cssample_set & cscontrol_set)
