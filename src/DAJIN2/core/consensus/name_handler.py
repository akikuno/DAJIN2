from __future__ import annotations

import re

from DAJIN2.core.consensus.consensus import ConsensusKey
from DAJIN2.utils.allele_handler import (
    AlleleComponents,
    build_allele_name,
    get_allele_base_name,
    parse_allele_name,
)


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
# Merge duplicate intact alleles
###########################################################


def merge_duplicate_intact_alleles(
    cons_sequences: dict[str, str],
    cons_percentages: dict[str, list],
    cons_midsv_tags: dict[str, list],
) -> tuple[dict[str, str], dict[str, list], dict[str, list]]:
    """
    Merge duplicate intact alleles that have identical sequences.

    Args:
        cons_sequences: Dictionary of allele names to consensus sequences
        cons_percentages: Dictionary of allele names to percentage data
        cons_midsv_tags: Dictionary of allele names to midsv tags

    Returns:
        Updated dictionaries with merged intact alleles
    """
    # Group alleles by their base name (e.g., "allele01_deletion01", "allele01_deletion_01_allele")
    allele_groups = {}
    for allele_name in cons_sequences:
        base_name = get_allele_base_name(allele_name)
        allele_groups.setdefault(base_name, []).append(allele_name)

    # Process each group to find and merge duplicate intact alleles
    merged_sequences = {}
    merged_percentages = {}
    merged_midsv_tags = {}

    for alleles in allele_groups.values():
        # Find all intact alleles in this group
        intact_alleles = [a for a in alleles if '_intact_' in a]

        if len(intact_alleles) > 1:
            # Check if all intact alleles have the same sequence
            sequences = [cons_sequences[a] for a in intact_alleles]
            if all(seq == sequences[0] for seq in sequences):
                # Merge the alleles
                # Calculate total percentage
                total_percent = sum(
                    float(parse_allele_name(a).percent) for a in intact_alleles
                )

                # Keep the first allele name but update percentage
                first_components = parse_allele_name(intact_alleles[0])
                merged_components = AlleleComponents(
                    allele_id=first_components.allele_id,
                    allele_name=first_components.allele_name,
                    suffix=first_components.suffix,
                    percent=f"{total_percent:.3f}"
                )
                merged_name = build_allele_name(merged_components)

                # Store merged data
                merged_sequences[merged_name] = cons_sequences[intact_alleles[0]]
                merged_percentages[merged_name] = cons_percentages[intact_alleles[0]]
                merged_midsv_tags[merged_name] = cons_midsv_tags[intact_alleles[0]]

                # Add non-intact alleles from this group
                for a in alleles:
                    if a not in intact_alleles:
                        merged_sequences[a] = cons_sequences[a]
                        merged_percentages[a] = cons_percentages[a]
                        merged_midsv_tags[a] = cons_midsv_tags[a]
            else:
                # Intact alleles have different sequences, keep them all
                for a in alleles:
                    merged_sequences[a] = cons_sequences[a]
                    merged_percentages[a] = cons_percentages[a]
                    merged_midsv_tags[a] = cons_midsv_tags[a]
        else:
            # No duplicate intact alleles, keep all alleles as is
            for a in alleles:
                merged_sequences[a] = cons_sequences[a]
                merged_percentages[a] = cons_percentages[a]
                merged_midsv_tags[a] = cons_midsv_tags[a]

    return merged_sequences, merged_percentages, merged_midsv_tags


###########################################################
# Add `NAME` key to RESULT_SAMPLE
###########################################################


def add_key_by_allele_name(clust_sample: list[dict], allele_names: dict[int, str]) -> list[dict]:
    for clust in clust_sample:
        label = clust["LABEL"]
        clust["NAME"] = allele_names[label]
    return clust_sample


def merge_result_sample_intact_alleles(result_sample: list[dict]) -> list[dict]:
    """
    Merge duplicate intact alleles in RESULT_SAMPLE.

    Args:
        result_sample: List of sample dictionaries with NAME and PERCENT fields

    Returns:
        Updated list with merged intact alleles
    """
    # Group by base allele name
    allele_groups = {}
    for item in result_sample:
        name = item["NAME"]
        base_name = get_allele_base_name(name)
        allele_groups.setdefault(base_name, []).append(item)

    # Process each group
    merged_results = []
    for items in allele_groups.values():
        # Find intact alleles
        intact_items = [item for item in items if '_intact_' in item["NAME"]]

        if len(intact_items) > 1:
            # Merge intact alleles
            total_percent = sum(item["PERCENT"] for item in intact_items)

            # Keep the first intact item and update its percentage
            merged_item = intact_items[0].copy()
            merged_item["PERCENT"] = total_percent

            # Update the name with the new percentage
            original_components = parse_allele_name(merged_item["NAME"])
            new_components = AlleleComponents(
                allele_id=original_components.allele_id,
                allele_name=original_components.allele_name,
                suffix=original_components.suffix,
                percent=f"{total_percent:.3f}"
            )
            merged_item["NAME"] = build_allele_name(new_components)

            merged_results.append(merged_item)

            # Add non-intact items
            for item in items:
                if item not in intact_items:
                    merged_results.append(item)
        else:
            # No duplicate intact alleles, keep all items
            merged_results.extend(items)

    return merged_results
