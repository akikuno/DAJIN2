def call_cssplit_diffloci(clust_sample, DIFFLOCI_ALLELES, MASKS_CONTROL):
    for cssplit in clust_sample:
        qname = cssplit["QNAME"]
        allele = cssplit["ALLELE"]
        sv = cssplit["SV"]
        cssplit = cssplit["CSSPLIT"].split(",")
        diffloci = DIFFLOCI_ALLELES[f"{allele}-{sv}"]
        maskloci = MASKS_CONTROL[allele]
