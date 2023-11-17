from __future__ import annotations
import re


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
    return ["=" + seq if seq.startswith("N") else seq for seq in cssplits]


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
