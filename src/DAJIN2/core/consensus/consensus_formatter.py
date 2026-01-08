from __future__ import annotations

import re
from collections import defaultdict
from typing import Any

from DAJIN2.utils.config import ConsensusKey


def scale_percentage(clust_sample_removed: list[dict]) -> list[dict]:
    total_percent = sum({c["LABEL"]: c["PERCENT"] for c in clust_sample_removed}.values())
    scaled_percent_by_label = {c["LABEL"]: round(c["PERCENT"] * 100 / total_percent, 3) for c in clust_sample_removed}
    for c in clust_sample_removed:
        c["PERCENT"] = scaled_percent_by_label[c["LABEL"]]
    return clust_sample_removed


def call_allele_name(
    cons_sequences: dict[ConsensusKey, str],
    FASTA_ALLELES: dict[str, str],
) -> tuple[dict[int, str], dict[int, str]]:
    """
    Nomenculature:
    allele{id}|{allele_name}|{allele_type}|{percent}%
    - allele01|control|intact|75%
    - allele02|unassigned|insertion|25%
    """
    digits = max(2, len(str(len(cons_sequences))))

    map_label_name = {}
    map_name_allele = {}
    for key, cons_seq in cons_sequences.items():
        allele_id = f"{key.label:0{digits}}"
        allele_name = key.allele
        if allele_name.endswith("_DAJIN2predicted"):
            match = re.match(r"(deletion|inversion|insertion)(\d+)", allele_name)
            allele_type = match.groups()[0]  # deletion, inversion, insertion
            allele_name = "unassigned"
        else:
            if cons_seq == FASTA_ALLELES[key.allele]:
                allele_type = "intact"
            else:
                allele_type = "indels"

        map_label_name[key.label] = f"allele{allele_id}|{allele_name}|{allele_type}|{key.percent}%"
        map_name_allele[f"allele{allele_id}|{allele_name}|{allele_type}|{key.percent}%"] = key.allele

    return map_label_name, map_name_allele


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
    clust_sample: list[dict], map_label_name: dict[int, str], label_before_to_after: dict[int, int]
) -> list[dict]:
    readnum_by_label = defaultdict(int)

    readnum_counted_label = set()
    for clust in clust_sample:
        old_label = clust["LABEL"]
        # Update LABEL, NAME, PERCENT
        new_label = label_before_to_after.get(old_label, old_label)
        new_name = map_label_name[new_label].replace("|", "_")
        new_percent = float(map_label_name[new_label].split("|")[-1].rstrip("%"))

        clust["LABEL"] = new_label
        clust["NAME"] = new_name
        clust["PERCENT"] = new_percent

        # Sum READNUM for merged labels
        if old_label not in readnum_counted_label:
            readnum_by_label[new_label] += clust["READNUM"]
            readnum_counted_label.add(old_label)

    # Finally, update READNUM
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


def _group_by_allele_and_seq(cons_sequences: dict[ConsensusKey, str]) -> dict[tuple, list[ConsensusKey]]:
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
def _merge_groups(groups: dict[tuple, list[ConsensusKey]]) -> tuple[list[dict], dict[int, int]]:
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
def _rank_labels_by_percent(merged_records: list[dict]) -> dict[int, int]:
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
def _compose_original_to_final_mapping(
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
    groups = _group_by_allele_and_seq(cons_sequences)
    merged_records, before_to_interim = _merge_groups(groups)
    interim_to_rank = _rank_labels_by_percent(merged_records)
    label_before_to_after = _compose_original_to_final_mapping(before_to_interim, interim_to_rank)

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
