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
# Convert cssplits to DNA sequence
###########################################################


def _revcomp_inversion(sequence: str) -> str:
    """Reverse complement only the lowercase portions (= inversion) of the DNA sequence."""
    # Define the complement dictionary
    complement = {"a": "t", "t": "a", "c": "g", "g": "c", "n": "n"}

    # Function to reverse complement a segment
    def reverse_complement(segment: str) -> str:
        return "".join(complement[nuc] for nuc in reversed(segment))

    # Use re to split sequence into lowercase segments and other parts
    parts = re.split(r"([a-z]+)", sequence)

    # Apply reverse complement to lowercase parts
    result = [reverse_complement(part) if part.islower() else part for part in parts]

    return "".join(result)


def convert_cssplits_to_sequence(cssplits: list[str]) -> str:
    sequence = []
    for tag in cssplits:
        if tag.startswith("-"):  # deletion
            pass
        elif tag.startswith("+"):  # insertion
            insertions, last_tag = tag.split("|")[:-1], tag.split("|")[-1]
            ins_seq = "".join([ins[1] for ins in insertions])
            sequence.append(ins_seq)
            if last_tag.startswith("-"):
                pass
            else:  # match or substitution
                last_seq = last_tag[-1]
            sequence.append(last_seq)
        else:  # match or substitution or inversion
            sequence.append(tag[-1])

    return _revcomp_inversion("".join(sequence))


###########################################################
# reverse complement to cssplits
###########################################################


def _reverse_cssplits(cssplits: list[str]) -> list[str]:
    for i, cs in enumerate(cssplits):
        if cs.startswith("+"):
            cssplits[i] = "+" + "|".join(cs.split("|")[::-1])
    return cssplits[::-1]


def _realign_insertion(cssplits: list[str]) -> list[str]:
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


def _complement_cssplit(cssplits: list[str]) -> list[str]:
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


def _add_match_operator_to_n(cssplits: list[str]) -> list[str]:
    """Add "=" (match operator) to the sequences with "N"."""
    cssplits_add_match_op = []
    for cs in cssplits:
        if cs.startswith("N") or cs.startswith("n"):
            cssplits_add_match_op.append("=" + cs)
        elif cs.startswith("+") and (cs[-1] == "N" or cs[-1] == "n"):
            cs_ins = cs.split("|")
            cs_last = "=" + cs_ins[-1]
            cs = "|".join(cs_ins[:-1]) + "|" + cs_last
            cssplits_add_match_op.append(cs)
        else:
            cssplits_add_match_op.append(cs)
    return cssplits_add_match_op


def _split_cssplits_by_delimiter(cssplits: list[str]) -> list[str]:
    cssplits_break = []
    for cs in cssplits:
        if cs.startswith("+"):
            cssplits_break.extend(cs.split("|"))
        else:
            cssplits_break.append(cs)
    return cssplits_break


def _combine_cssplits_by_prefix(cssplits: list[str]) -> list[str]:
    if not cssplits:
        return []

    combined_cssplits = []
    prev_prefix: str = cssplits[0][0]
    current_cs: list[str] = [cssplits[0][1:]]

    for item in cssplits[1:]:
        current_prefix, cs = item[0], item[1:]
        if prev_prefix == current_prefix:
            if current_prefix == "*":
                current_cs.append(item)
            else:
                current_cs.append(cs)
        else:
            combined_cssplits.append(prev_prefix + "".join(current_cs))
            prev_prefix = current_prefix
            current_cs = [cs]

    combined_cssplits.append(prev_prefix + "".join(current_cs))
    return combined_cssplits


def _standardize_case(cssplits: list[str]) -> list[str]:
    """Standardize the case of characters based on mutation types."""
    transformations = {"*": str.lower, "-": str.lower, "+": str.lower, "~": str.lower, "=": str.upper}

    transformed = []
    for cs in cssplits:
        prefix = cs[0]
        trans_func = transformations.get(prefix)
        transformed.append(prefix + trans_func(cs[1:]))

    return transformed


def convert_cssplits_to_cstag(cssplits: list[str]) -> str:
    cssplits_matched_op = _add_match_operator_to_n(cssplits)
    cssplits_splitted = _split_cssplits_by_delimiter(cssplits_matched_op)
    cssplits_combined = _combine_cssplits_by_prefix(cssplits_splitted)
    return "".join(_standardize_case(cssplits_combined))


###########################################################
# convert cssplits to sequence
###########################################################


def call_sequence(cons_percentage: list[dict[str, float]]) -> str:
    """convert position weight matrix (cons_pergentage) to sequence"""

    consensus_sequence = []
    for cons_per in cons_percentage:
        cssplits = max(cons_per, key=cons_per.get)
        cs_tag = convert_cssplits_to_cstag([cssplits])
        seq = cstag.to_sequence(cs_tag)
        consensus_sequence.append(seq)
    return "".join(consensus_sequence)
