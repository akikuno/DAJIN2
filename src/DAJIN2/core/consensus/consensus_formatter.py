from __future__ import annotations

import re
from collections import defaultdict
from dataclasses import dataclass
from typing import Any


@dataclass(frozen=True)
class ConsensusKey:
    allele: str
    label: int
    percent: float


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
    sv_threshold: int = 50,
) -> dict[int, str]:
    digits = max(2, len(str(len(cons_percentages))))
    exists_sv = [detect_sv(cons_per, sv_threshold) for cons_per in cons_percentages.values()]

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


def update_key_by_allele_name(cons: dict[ConsensusKey, Any], allele_names: dict[int, str]) -> dict:
    cons_update = {}
    for consensuskey, value in cons.items():
        allele_name = allele_names[consensuskey.label]
        cons_update[allele_name] = value
    return cons_update


###########################################################
# Update `LABEL`, `PERCENT`, and `NAME` key to RESULT_SAMPLE
# After clustering, alleles with the same consensus sequence were merged.
# Therefore, Label and Percent have changed, and they are updated here.
###########################################################


def update_label_percent_readnum_name(
    clust_sample: list[dict], allele_names: dict[int, str], label_before_to_after: dict[int, int]
) -> list[dict]:
    readnum_by_label = defaultdict(int)

    readnum_counted_label = set()
    for clust in clust_sample:
        old_label = clust["LABEL"]
        new_label = label_before_to_after.get(old_label, old_label)
        new_name = allele_names[new_label]
        new_percent = float(new_name.split("_")[-1].rstrip("%"))

        clust["LABEL"] = new_label
        clust["NAME"] = new_name
        clust["PERCENT"] = new_percent

        if old_label not in readnum_counted_label:
            readnum_by_label[new_label] += clust["READNUM"]
            readnum_counted_label.add(old_label)

    # Overwrite READNUM in the second loop
    for clust in clust_sample:
        clust["READNUM"] = readnum_by_label[clust["LABEL"]]

    return clust_sample


###########################################################
# Merge identical alleles when the same consensus sequence
# has been assigned different labels.
# While consolidating, sum percentages and reassign labels
# as 1, 2, 3, ... in descending order of total percentage.
###########################################################

# -----------------------------
# 1) Grouping
# -----------------------------


def group_by_allele_and_seq(cons_sequences: dict[ConsensusKey, str]) -> dict[tuple, list[ConsensusKey]]:
    """
    Group ConsensusKey objects by (allele, sequence).

    Returns:
        {(allele, seq): [ConsensusKey, ...], ...}
    """
    groups: dict[tuple, list[ConsensusKey]] = defaultdict(list)
    for key, seq in cons_sequences.items():
        groups[(key.allele, seq)].append(key)
    return dict(groups)


# -----------------------------
# 2) Merge (pick representative & sum percentages & build before→rep mapping)
# -----------------------------
def merge_groups(groups: dict[tuple, list[ConsensusKey]]) -> tuple[list[dict], dict[int, int]]:
    """
    For each group:
        - Create an interim_label (a temporary label after merging but before sorting by percentage).
        - Sum the percentages.
        - Build before_to_interim[label_before] = interim_label (exclude the representative itself).

    Returns:
        merged_records: [{"allele": str, "seq": str, "interim_label": int, "percent": float}, ...]
        before_to_interim: {label_before: interim_label, ...}
    """
    merged_records: list[dict] = []
    before_to_interim: dict[int, int] = {}

    for (allele, seq), keys in groups.items():
        interim_label = keys[0].label  # choose any key as the representative (deterministic: first seen)
        total_percent = sum(float(k.percent) for k in keys)
        for k in keys:
            if k.label != interim_label:
                before_to_interim[k.label] = interim_label
        merged_records.append({"allele": allele, "seq": seq, "interim_label": interim_label, "percent": total_percent})
    return merged_records, before_to_interim


# -----------------------------
# 3) Ranking (representative label → scores 1..N)
# -----------------------------
def rank_labels_by_percent(merged_records: list[dict]) -> dict[int, int]:
    """
    Reassign labels 1, 2, 3, ... in descending order of the summed percentage.

    Returns:
        interim_to_rank: {interim_label: rank_label, ...}
    """
    ordered = sorted(merged_records, key=lambda r: (-float(r["percent"]), int(r["interim_label"])))
    interim_to_rank: dict[int, int] = {r["interim_label"]: i + 1 for i, r in enumerate(ordered)}
    return interim_to_rank


# -----------------------------
# 4) Compose mapping (original label → final label)
# -----------------------------
def compose_original_to_final_mapping(
    before_to_interim: dict[int, int], interim_to_rank: dict[int, int]
) -> dict[int, int]:
    """
    Build a mapping from original labels (including representatives) to final labels (ranks).
    - Map each representative via interim_to_rank to its rank.
    - Map each non-representative via before_to_interim → interim_to_rank to its rank.

    Returns:
        original_to_final: {original_label: final_rank_label, ...}
    """
    original_to_final: dict[int, int] = {}

    # Representatives: map directly to their rank
    for interim_label, rank in interim_to_rank.items():
        original_to_final[interim_label] = rank

    # Non-representatives: map via representative → rank
    for before_label, interim_label in before_to_interim.items():
        original_to_final[before_label] = interim_to_rank[interim_label]

    return original_to_final


# -----------------------------
# Convenience wrapper: composed pipeline
# -----------------------------
def merge_duplicated_cons_sequences(
    cons_sequences: dict[ConsensusKey, str],
) -> tuple[dict[ConsensusKey, str], dict[int, int]]:
    """
    1) Group → 2) Merge → 3) Rank → 4) Compose mapping → 5) Return dict with relabeled keys.

    Returns:
        cons_sequences_merged:
            {ConsensusKey(allele, label=rank, percent=merged_percent): seq, ...}
        label_before_to_after:
            {original_label: final_rank_label, ...}
    """
    groups = group_by_allele_and_seq(cons_sequences)
    merged_records, before_to_interim = merge_groups(groups)
    interim_to_rank = rank_labels_by_percent(merged_records)
    label_before_to_after = compose_original_to_final_mapping(before_to_interim, interim_to_rank)

    # Rebuild the dict with rank labels
    cons_sequences_merged: dict[ConsensusKey, str] = {}
    for rec in merged_records:
        rank_label = interim_to_rank[rec["interim_label"]]
        key = ConsensusKey(allele=rec["allele"], label=rank_label, percent=rec["percent"])
        cons_sequences_merged[key] = rec["seq"]

    return cons_sequences_merged, label_before_to_after


def merge_duplicated_cons_others(
    cons_others: dict[ConsensusKey, Any],
    cons_sequences_merged: dict[ConsensusKey, str],
    label_before_to_after: dict[int, int],
) -> dict[ConsensusKey, Any]:
    unique_alleles = {(item.allele, item.label): item.percent for item in cons_sequences_merged}

    cons_others_merged = {}
    for key, value in cons_others.items():
        new_label = label_before_to_after[key.label]
        if (key.allele, new_label) not in unique_alleles:
            continue
        percent = unique_alleles[(key.allele, new_label)]
        new_key = ConsensusKey(allele=key.allele, label=new_label, percent=percent)
        cons_others_merged[new_key] = value

    return cons_others_merged
