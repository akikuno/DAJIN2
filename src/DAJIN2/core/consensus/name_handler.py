from __future__ import annotations

import re
from pathlib import Path

from DAJIN2.core.consensus.consensus import ConsensusKey


def scale_percentage(clust_sample_removed: list[dict]) -> list[dict]:
    total_percent = sum({c["LABEL"]: c["PERCENT"] for c in clust_sample_removed}.values())
    scaled_percent_by_label = {c["LABEL"]: round(c["PERCENT"] * 100 / total_percent, 3) for c in clust_sample_removed}
    for c in clust_sample_removed:
        c["PERCENT"] = scaled_percent_by_label[c["LABEL"]]
    return clust_sample_removed


def detect_sv(cons_per: list[dict[str, float]], threshold: int = 50) -> bool:
    cons_midsv_tag: str = "".join([max(tag, key=tag.get) for tag in cons_per])

    # TODO: "N" should be replaced with "=N" in the next MIDSV update
    patterns = [
        rf"N{{{threshold}}}",  # Consecutive "N" exceeding the threshold
        rf"(\+[ACGTN]\|){{{threshold}}}",  # Insertions
        rf"(\-[ACGTN]){{{threshold}}}",  # Deletions
        rf"(\*[ACGTN][ACGTN]){{{threshold}}}",  # Substitutions
        r"[acgtn]",  # Inversions (lowercase nucleotides)
    ]

    return any(re.search(pattern, cons_midsv_tag) for pattern in patterns)


def format_allele_label(label: int, total_labels: int) -> str:
    digits = max(2, len(str(total_labels)))  # minimum of 2 digits (01, 02, 03...)
    return f"{label:0{digits}}"


def determine_suffix(cons_seq: str, fasta_allele: str, is_sv: bool) -> str:
    if is_sv:
        return "SV"
    elif cons_seq == fasta_allele:
        return "intact"
    else:
        return "indels"


def generate_allele_mapping(alleles: list[str]) -> dict[str, str]:
    # Define the mapping and groups
    groups = {}
    # Group alleles by prefix (deletion, inversion, insertion)
    for allele in alleles:
        match = re.match(r"(deletion|inversion|insertion)(\d+)", allele)
        if match:
            prefix, _ = match.groups()
            if prefix not in groups:
                groups[prefix] = []
            groups[prefix].append(allele)

    # Sort each group by percent (descending) and assign new numbers
    allele_mapping = {}
    for prefix, group in groups.items():
        digits = max(2, len(str(len(group))))
        allele_id = 1
        for allele in group:
            if allele not in allele_mapping:
                new_allele = f"{prefix}{(allele_id):0{digits}}"
                allele_mapping[allele] = new_allele
                allele_id += 1

    return allele_mapping


def call_allele_name(
    tempdir: Path | str,
    sample_name: str,
    cons_sequences: dict[ConsensusKey, str],
    cons_percentages: dict[ConsensusKey, list],
    FASTA_ALLELES: dict[str, str],
    threshold: int = 50,
) -> dict[int, str]:
    digits = max(2, len(str(len(cons_percentages))))
    exists_sv = [detect_sv(cons_per, threshold) for cons_per in cons_percentages.values()]

    sorted_keys = sorted(cons_percentages, key=lambda x: x.percent, reverse=True)
    alleles = [key.allele for key in sorted_keys]
    allele_mapping = generate_allele_mapping(alleles)

    allele_names = {}
    for is_sv, (keys, cons_seq) in zip(exists_sv, cons_sequences.items()):
        allele_id = f"{keys.label:0{digits}}"
        allele_name = allele_mapping.get(keys.allele, keys.allele)
        suffix = determine_suffix(cons_seq, FASTA_ALLELES[keys.allele], is_sv)
        allele_names[keys.label] = f"allele{allele_id}_{allele_name}_{suffix}_{keys.percent}%"

    return allele_names


###########################################################
# Replace new allele names to the consensus dictionary
###########################################################


def update_key_by_allele_name(cons: dict, allele_names: dict[int, str]) -> dict:
    cons_update = {}
    for key in cons:
        old_allele = cons[key]
        new_allele = allele_names[key.label]
        cons_update[new_allele] = old_allele
    return cons_update


###########################################################
# Add `NAME` key to RESULT_SAMPLE
###########################################################


def add_key_by_allele_name(clust_sample: list[dict], allele_names: dict[int, str]) -> list[dict]:
    for clust in clust_sample:
        label = clust["LABEL"]
        clust["NAME"] = allele_names[label]
    return clust_sample
