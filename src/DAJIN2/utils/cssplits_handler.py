from __future__ import annotations
import re
import cstag


def find_n_boundaries(cssplits: list[str]) -> tuple[int, int]:
    """Find the boundaries of contiguous Ns which aren't at the ends."""

    # Find the left boundary
    left_idx_n = 0
    for char in cssplits:
        if char != "N":
            break
        left_idx_n += 1

    # Find the right boundary
    right_idx_n = len(cssplits) - 1
    for char in reversed(cssplits):
        if char != "N":
            break
        right_idx_n -= 1

    return left_idx_n - 1, right_idx_n + 1


###########################################################
# reverse complement to cssplits
###########################################################


def _reverse_cssplits(cssplits: list) -> list:
    for i, cs in enumerate(cssplits):
        if cs.startswith("+"):
            cssplits[i] = "+" + "|".join(cs.split("|")[::-1])
    return cssplits[::-1]


def _realign_insertion(cssplits: list) -> list:
    for i, cs in enumerate(cssplits):
        if not cs.startswith("+"):
            continue
        if re.search(rf"[{cs[1]}]", "[ACGTacgt]"):
            continue
        if i + 1 == len(cssplits):
            continue
        cs_current = cs.split("|")
        cssplits[i] = cs_current[0].replace("+", "")
        cssplits[i + 1] = "|".join(c[0] + c[-1] for c in cs_current[1:]) + "|" + cssplits[i + 1]
    return cssplits


def _complement_cssplit(cssplits: list) -> list:
    comp = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N", "a": "t", "c": "g", "g": "c", "t": "a", "n": "n"}
    for i, cs in enumerate(cssplits):
        op = cs[0]
        if op == "*":
            cssplits[i] = op + comp[cs[1]] + comp[cs[2]]
        elif op == "+":
            cssplits[i] = "|".join(c[0] + comp[c[-1]] for c in cs.split("|"))
        else:  # Match or Deletion or N
            cssplits[i] = op + comp[cs[-1]]
        cssplits[i] = cssplits[i].replace("NN", "N")
        cssplits[i] = cssplits[i].replace("nn", "n")
    return cssplits


def revcomp_cssplits(cssplits: list[str]) -> list[str]:
    cssplits_reversed = _reverse_cssplits(cssplits)
    cssplits_realigned = _realign_insertion(cssplits_reversed)
    cssplits_revcomped = _complement_cssplit(cssplits_realigned)
    return cssplits_revcomped


###########################################################
# convert cssplits to cstag
###########################################################


def add_match_operator_to_n(cssplits: list[str]) -> list[str]:
    """Add "=" (match operator) to the sequences that start with "N"."""
    return ["=" + seq if seq.startswith("N") or seq.startswith("n") else seq for seq in cssplits]


def format_insertion(cs: str) -> str:
    """Reformat insertion sequence by consolidating variants."""
    # Remove "+" and split by "|"
    variants = cs.replace("+", "").split("|")
    # Extract the last character from all but the last variant
    consolidated_variants = "".join([variant[-1] for variant in variants[:-1]])
    # Combine consolidated variants with the last variant
    formatted_insertion = "+" + consolidated_variants + variants[-1]

    return formatted_insertion


def concatenate_cssplits(cssplits: list[str]) -> str:
    """Concatenate list of sequences based on certain variants."""
    if not cssplits:
        return ""
    concatenated = []
    prev = cssplits[0]
    if prev[0] == "+":
        concatenated.append(format_insertion(prev))
    else:
        concatenated.append(prev)
    for prev, current in zip(cssplits, cssplits[1:]):
        if prev[0] == current[0] and current[0] in {"=", "-"}:
            concatenated.append(current[1:])
        elif current[0] == "+":
            concatenated.append(format_insertion(current))
        else:
            concatenated.append(current)

    return "".join(concatenated)


def standardize_case(sequence: str) -> str:
    """Standardize the case of characters based on preceding variants."""
    transformations = {"*": str.lower, "-": str.lower, "+": str.lower, "~": str.lower, "=": str.upper}

    transformed = []
    i = 0
    while i < len(sequence):
        char = sequence[i]
        trans_func = transformations.get(char)
        if trans_func:
            start = i + 1
            while start < len(sequence) and sequence[start].isalpha():
                start += 1
            transformed.append(char + trans_func(sequence[i + 1 : start]))
            i = start - 1
        else:
            transformed.append(char)
        i += 1

    return "".join(transformed)


def convert_cssplits_to_cstag(cssplits: list[str]) -> str:
    cssplits = add_match_operator_to_n(cssplits)
    cssplits_concatenated = concatenate_cssplits(cssplits)
    return standardize_case(cssplits_concatenated)


# convert position weight matrix (cons_pergentage) to sequence


def call_sequence(cons_percentage: list[dict[str, float]], sep: str = "") -> str:
    consensus_sequence = []
    for cons_per in cons_percentage:
        cssplits = max(cons_per, key=cons_per.get)
        cs_tag = convert_cssplits_to_cstag([cssplits])
        seq = cstag.to_sequence(cs_tag)
        consensus_sequence.append(seq)
    return f"{sep}".join(consensus_sequence)


###########################################################
# detect_insertion_within_deletion
###########################################################


def is_start_of_deletion(i: int, cssplits: list[str], start_del_count: int = 3) -> bool:
    """
    Determine if the current index is the start of a deletion.
    If there are consecutive deletions from the beginning, it is determined that there is a cluster of deletions.
    """
    if i + start_del_count > len(cssplits):
        return False
    return all([cs.startswith("-") for cs in cssplits[i : i + start_del_count]])


def is_within_deletion(i: int, cssplits: list[str], distance: int = 10) -> bool:
    """Count matches that are continuous for `distance` bases, and if there is a deletion in between, it is determined to be inside a cluster of deletions."""
    exist_deletion = False
    total_num_match = 0
    j = 1
    while True:
        if i + j >= len(cssplits) or total_num_match > distance:
            exist_deletion = False
            break
        current = cssplits[i + j]
        if current.startswith("="):
            total_num_match += 1
        else:
            total_num_match = 0
            if current.startswith("-"):
                exist_deletion = True
                break
        j += 1
    return exist_deletion


def adjust_cs_insertion(cs: str) -> str:
    if cs.startswith("-"):
        return None
    elif cs.startswith("+"):
        ins = "|".join(cs.split("|")[:-1])
        end_cs = cs.split("|")[-1][-1]
        cs_insertion = f"{ins}|+{end_cs}|"
    else:
        cs_insertion = f"+{cs[-1]}|"
    return cs_insertion


def detect_insertion_within_deletion(cssplits: str, start_del_count: int = 3, distance: int = 10) -> str:
    cssplits = cssplits.split(",")
    cssplits_updated = cssplits.copy()

    insertion = []
    within_deletion = False
    for i, cs in enumerate(cssplits):
        if cs.startswith("-"):
            if within_deletion is False and is_start_of_deletion(i, cssplits, start_del_count) is False:
                continue
            within_deletion = is_within_deletion(i, cssplits, distance)
            if within_deletion is False:
                cssplits_updated[i + 1] = "".join(insertion) + cssplits_updated[i + 1]
                insertion = []
        if within_deletion is False:
            continue
        if not cs.startswith("-"):
            if cs.startswith("*"):
                cssplits_updated[i] = f"-{cs[1]}"
            elif cs.startswith("+") and cs.split("|")[-1].startswith("*"):
                cssplits_updated[i] = f'-{cs.split("|")[-1][1]}'
            else:
                cssplits_updated[i] = f"-{cs[-1]}"
            adjusted_cs = adjust_cs_insertion(cs)
            insertion.append(adjusted_cs)

    return ",".join(cssplits_updated)
