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


def is_start_of_consecutive_deletions(cssplits: list[str], i: int, del_range: int = 3) -> bool:
    """Determine if the current index is the start of a consective deletion."""
    if i + del_range > len(cssplits):
        return False

    for j in range(del_range):
        if not cssplits[i + j].startswith("-"):
            return False

    return True


def get_range_of_inclusive_deletions(cssplits: list[str], del_range: int = 3, distance: int = 10) -> set[int]:
    inclusive_deletion_range: set[int] = set()
    n: int = len(cssplits)
    i: int = 0
    while i < n:
        # Check for start of deletion
        if cssplits[i].startswith("-"):
            start: int = i
            while i < n and cssplits[i].startswith("-"):
                i += 1
            # We have found the end of a deletion range
            end: int = i
            # Check if deletion length is enough to be considered

            if end - start >= del_range:
                # Now check for matching sequence within 'distance'
                match_count: int = 0
                while i < n and match_count <= distance:
                    if cssplits[i].startswith("="):
                        match_count += 1
                    elif cssplits[i].startswith("-"):
                        # If another deletion is found, reset match count and update range
                        next_start: int = i
                        while i < n and cssplits[i].startswith("-"):
                            i += 1
                            next_end: int = i
                        if next_end - next_start >= del_range:
                            inclusive_deletion_range.update(range(start, next_end))

                            start = next_start
                            end = next_end
                            match_count = 0
                            break
                    i += 1

        i += 1

    return inclusive_deletion_range


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


def reallocate_insertion_within_deletion(cssplit: str, del_range: int = 3, distance: int = 10) -> str:
    """
    Since the mapping in minimap2 is local alignment, insertion bases within large deletions may be partially mapped to the reference genome and not detected as insertion bases. Therefore, update cssplits to detect insertions within large deletions as insertions.
    """
    cssplits = cssplit.split(",")
    cssplits_updated = cssplits.copy()

    range_of_inclusive_deletion: set[int] = get_range_of_inclusive_deletions(cssplits, del_range, distance)

    insertion_within_deletion = []
    for i, cs in enumerate(cssplits):
        if i in range_of_inclusive_deletion:
            if cs.startswith("-"):
                continue
            if cs.startswith("N") or cs.startswith("n"):
                cs_new = cs
            elif cs.startswith("*"):
                cs_new = f"-{cs[1]}"
            elif cs.startswith("+") and cs.split("|")[-1].startswith("*"):
                cs_new = f'-{cs.split("|")[-1][1]}'
            else:
                cs_new = f"-{cs[-1]}"
            cssplits_updated[i] = cs_new

            insertion_within_deletion.append(adjust_cs_insertion(cs))
        else:
            if not insertion_within_deletion:
                continue
            cssplits_updated[i] = "".join(insertion_within_deletion) + cs
            insertion_within_deletion = []

    return ",".join(cssplits_updated)
