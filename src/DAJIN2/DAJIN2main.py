from __future__ import annotations

import warnings

warnings.simplefilter("ignore")
import hashlib
from collections import defaultdict
from copy import deepcopy
from itertools import groupby
from pathlib import Path

import midsv

from src.DAJIN2 import classification, clustering
from src.DAJIN2.consensus import module_consensus as consensus
from src.DAJIN2.preprocess import check_inputs, format_inputs, mappy_align


def assignment(arguments: dict, key: str, value: str):
    try:
        return arguments[key]
    except KeyError:
        return value


def main(arguments: dict) -> None:
    SAMPLE = arguments["sample"]
    CONTROL = arguments["control"]
    ALLELE = arguments["allele"]
    OUTPUT = assignment(arguments, "output", "DAJIN_results")
    GENOME = assignment(arguments, "genome", "")
    THREADS = int(assignment(arguments, "threads", "1"))
    ##########################################################
    # Check inputs
    ##########################################################
    check_inputs.check_files(SAMPLE, CONTROL, ALLELE)
    IS_CACHE_CONTROL = check_inputs.is_cache_control(CONTROL, OUTPUT)
    IS_CACHE_GENOME = check_inputs.is_cache_genome(GENOME, OUTPUT, IS_CACHE_CONTROL)
    if GENOME and not IS_CACHE_GENOME:
        UCSC_URL, GOLDENPATH_URL = check_inputs.check_and_fetch_genome(GENOME)
    ##########################################################
    # Format inputs
    ##########################################################
    SAMPLE_NAME = format_inputs.extract_basename(SAMPLE)
    CONTROL_NAME = format_inputs.extract_basename(CONTROL)
    DICT_ALLELE = format_inputs.dictionize_allele(ALLELE)

    if GENOME:
        path_genome_coodinates = Path(OUTPUT, ".tempdir", "cache", "genome_coodinates.jsonl")
        path_chrome_size = Path(OUTPUT, ".tempdir", "cache", "chrome_size.txt")
        if IS_CACHE_GENOME:
            GENOME_COODINATES = midsv.read_jsonl(path_genome_coodinates)[0]
            CHROME_SIZE = int(path_chrome_size.read_text())
        else:
            GENOME_COODINATES = format_inputs.fetch_coodinate(GENOME, UCSC_URL, DICT_ALLELE["control"])
            CHROME_SIZE = format_inputs.fetch_chrom_size(GENOME_COODINATES["chr"], GENOME, GOLDENPATH_URL)
            # Save info to the cache directory
            Path(OUTPUT, ".tempdir", "cache", "genome_symbol.txt").write_text(GENOME)
            midsv.write_jsonl([GENOME_COODINATES], path_genome_coodinates)
            path_chrome_size.write_text(str(CHROME_SIZE))

    format_inputs.make_directories(OUTPUT)

    ################################################################################
    # Export fasta files as single-FASTA format
    ################################################################################
    # TODO: use yeild, not export
    for identifier, sequence in DICT_ALLELE.items():
        contents = "\n".join([">" + identifier, sequence]) + "\n"
        output_fasta = Path(OUTPUT, ".tempdir", "fasta", f"{identifier}.fasta")
        output_fasta.write_text(contents)

    ###############################################################################
    # Mapping with minimap2/mappy
    ###############################################################################
    for input_fasta in Path(OUTPUT, ".tempdir", "fasta").glob("*.fasta"):
        identifier = input_fasta.stem
        if identifier not in set(DICT_ALLELE.keys()):
            continue
        if not IS_CACHE_CONTROL:
            fastq, fastq_name = CONTROL, CONTROL_NAME
            sam = mappy_align.to_sam(str(input_fasta), fastq)
            output_sam = Path(OUTPUT, ".tempdir", "sam", f"{fastq_name}_{identifier}.sam")
            output_sam.write_text("\n".join(sam))
        fastq, fastq_name = SAMPLE, SAMPLE_NAME
        sam = mappy_align.to_sam(str(input_fasta), fastq)
        output_sam = Path(OUTPUT, ".tempdir", "sam", f"{fastq_name}_{identifier}.sam")
        output_sam.write_text("\n".join(sam))

    ########################################################################
    # MIDSV conversion
    ########################################################################
    if not IS_CACHE_CONTROL:
        for sampath in Path(OUTPUT, ".tempdir", "sam").glob(f"{CONTROL_NAME}*"):
            identifier = sampath.stem.split("_")[1]
            if identifier not in set(DICT_ALLELE.keys()):
                continue
            sam = midsv.read_sam(sampath)
            midsv_jsonl = midsv.transform(sam, midsv=False, cssplit=True, qscore=False)
            output_jsonl = Path(OUTPUT, ".tempdir", "midsv", f"{sampath.stem}.jsonl")
            midsv.write_jsonl(midsv_jsonl, output_jsonl)

    for sampath in Path(OUTPUT, ".tempdir", "sam").glob(f"{SAMPLE_NAME}*"):
        identifier = sampath.stem.split("_")[1]
        if identifier not in set(DICT_ALLELE.keys()):
            continue
        sam = midsv.read_sam(sampath)
        midsv_jsonl = midsv.transform(sam, midsv=False, cssplit=True, qscore=False)
        output_jsonl = Path(OUTPUT, ".tempdir", "midsv", f"{sampath.stem}.jsonl")
        midsv.write_jsonl(midsv_jsonl, output_jsonl)

    ###############################################################################
    # Cashe inputs (control)
    ###############################################################################

    if not IS_CACHE_CONTROL:
        control_hash = Path(CONTROL).read_bytes()
        control_hash = hashlib.sha256(control_hash).hexdigest()
        PATH_CACHE_HASH = Path(OUTPUT, ".tempdir", "cache", "control_hash.txt")
        PATH_CACHE_HASH.write_text(str(control_hash))

    ########################################################################
    # Classify alleles
    ########################################################################
    path_midsv = Path(OUTPUT, ".tempdir", "midsv").glob(f"{SAMPLE_NAME}*")
    classif_sample = classification.classify_alleles(path_midsv, SAMPLE_NAME)
    path_midsv = Path(OUTPUT, ".tempdir", "midsv").glob(f"{CONTROL_NAME}*")
    classif_control = classification.classify_alleles(path_midsv, CONTROL_NAME)
    ########################################################################
    # Detect Structural variants
    ########################################################################
    for classifs in [classif_sample, classif_control]:
        for classif in classifs:
            classif["SV"] = classification.detect_sv(classif["CSSPLIT"], threshold=50)

    ########################################################################
    # Clustering
    ########################################################################
    # -----------------------------------------------------------------------
    # Extract significantly different base loci between Sample and Control
    # -----------------------------------------------------------------------
    dict_cssplit_control = defaultdict(list[dict])
    for ALLELE in DICT_ALLELE.keys():
        path_control = Path(OUTPUT, ".tempdir", "midsv", f"{CONTROL_NAME}_{ALLELE}.jsonl")
        cssplit_control = [cs["CSSPLIT"] for cs in midsv.read_jsonl(path_control)]
        dict_cssplit_control[ALLELE] = cssplit_control

    classif_sample.sort(key=lambda x: (x["ALLELE"], x["SV"]))
    diffloci_by_alleles = defaultdict(list[dict])
    for (ALLELE, SV), group in groupby(classif_sample, key=lambda x: (x["ALLELE"], x["SV"])):
        cssplit_sample = [record["CSSPLIT"] for record in group]
        cssplit_control = dict_cssplit_control[ALLELE]
        sequence = DICT_ALLELE[ALLELE]
        diffloci = clustering.screen_different_loci(
            cssplit_sample, cssplit_control, sequence, alpha=0.01, threshold=0.05
        )
        diffloci_by_alleles[f'{{"ALLELE": "{ALLELE}", "SV": {SV}}}'] = diffloci

    # -----------------------------------------------------------------------
    # Clustering
    # -----------------------------------------------------------------------
    labels = []
    label_start = 1
    for (ALLELE, SV), group in groupby(classif_sample, key=lambda x: (x["ALLELE"], x["SV"])):
        key = f'{{"ALLELE": "{ALLELE}", "SV": {SV}}}'
        cssplit_sample = [g["CSSPLIT"] for g in group]
        diffloci = diffloci_by_alleles[key]
        scores = list(clustering.make_scores(cssplit_sample, diffloci))
        if any(scores):
            labels += [label + label_start for label in clustering.clustering(scores).tolist()]
        else:
            labels += [label_start] * len(cssplit_sample)
        label_start = len(set(labels)) + 1

    clust_sample = deepcopy(classif_sample)
    for clust, label in zip(clust_sample, labels):
        clust["LABEL"] = label
        del clust["CSSPLIT"]

    n_sample = len(clust_sample)
    d = defaultdict(int)
    for cs in clust_sample:
        d[cs["LABEL"]] += 1 / n_sample

    d_per = {key: round(val * 100, 1) for key, val in d.items()}

    for cs in clust_sample:
        cs["PERCENT"] = d_per[cs["LABEL"]]

    # Allocate new labels by PERCENT
    clust_sample.sort(key=lambda x: (-x["PERCENT"], x["LABEL"]))
    new_label = 1
    prev_label = clust_sample[0]["LABEL"]
    for cs in clust_sample:
        if prev_label != cs["LABEL"]:
            new_label += 1
        prev_label = cs["LABEL"]
        cs["LABEL"] = new_label

    ########################################################################
    # Consensus call
    ########################################################################
    path = Path(OUTPUT, ".tempdir", "midsv", f"{CONTROL_NAME}_control.jsonl")
    cssplit_control = midsv.read_jsonl(path)
    path = Path(OUTPUT, ".tempdir", "midsv", f"{SAMPLE_NAME}_control.jsonl")
    cssplit_sample = midsv.read_jsonl(path)
    cssplit_sample = consensus.join_listdicts(clust_sample, cssplit_sample, key="QNAME")
    cssplit_sample.sort(key=lambda x: (x["ALLELE"], x["SV"], x["LABEL"]))
    cons_percentage = defaultdict(list)
    cons_sequence = defaultdict(list)
    for keys, cssplits in groupby(cssplit_sample, key=lambda x: (x["ALLELE"], x["SV"], x["LABEL"])):
        cssplits = list(cssplits)
        cons_per = consensus.call_percentage(cssplits, cssplit_control)
        cons_seq = consensus.call_sequence(cons_per)
        allele_name = consensus.call_allele_name(keys, cons_seq, DICT_ALLELE)
        cons_percentage[allele_name] = cons_per
        cons_sequence[allele_name] = cons_seq
        for cs in cssplits:
            cs["NAME"] = allele_name

    # ----------------------------------------------------------
    # Conseusns Report：FASTA/HTML/VCF
    # ----------------------------------------------------------
    # FASTA
    for header, cons_seq in cons_sequence.items():
        cons_fasta = consensus.to_fasta(header, cons_seq)
        Path(f"{OUTPUT}/.tempdir/reports/{SAMPLE_NAME}_{header}.fasta").write_text(cons_fasta)

    # HTML
    for header, cons_per in cons_percentage.items():
        cons_html = consensus.to_html(header, cons_per)
        Path(f"{OUTPUT}/.tempdir/reports/{SAMPLE_NAME}_{header}.html").write_text(cons_html)

    # VCF
    # working in progress

    RESULT_SAMPLE = deepcopy(cssplit_sample)
    for res in RESULT_SAMPLE:
        del res["RNAME"]
        del res["CSSPLIT"]
    return RESULT_SAMPLE
