from __future__ import annotations

import re

import cstag


def is_n_tag(tag: str) -> bool:
    return tag.endswith("N") or tag.endswith("n")


def find_n_boundaries(midsv_tags: list[str]) -> tuple[int, int]:
    """Find the boundaries of contiguous Ns which aren't at the ends."""

    # Find the left boundary
    left_idx_n = 0
    for char in midsv_tags:
        if not is_n_tag(char):
            break
        left_idx_n += 1

    # Find the right boundary
    right_idx_n = len(midsv_tags) - 1
    for char in reversed(midsv_tags):
        if not is_n_tag(char):
            break
        right_idx_n -= 1

    return left_idx_n - 1, right_idx_n + 1


###########################################################
# Convert midsv_tags to DNA sequence
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


def convert_midsvs_to_sequence(midsv_tags: list[str]) -> str:
    sequence = []
    for tag in midsv_tags:
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
# reverse complement to midsv_tags
###########################################################


def _reverse_midsvs(midsv_tags: list[str]) -> list[str]:
    for i, cs in enumerate(midsv_tags):
        if cs.startswith("+"):
            midsv_tags[i] = "+" + "|".join(cs.split("|")[::-1])
    return midsv_tags[::-1]


def _realign_insertion(midsv_tags: list[str]) -> list[str]:
    for i, cs in enumerate(midsv_tags):
        if not cs.startswith("+"):
            continue
        if re.search(rf"[{cs[1]}]", "[ACGTacgt]"):
            continue
        if i + 1 == len(midsv_tags):
            continue
        cs_current = cs.split("|")
        midsv_tags[i] = cs_current[0].replace("+", "")
        midsv_tags[i + 1] = "|".join(c[0] + c[-1] for c in cs_current[1:]) + "|" + midsv_tags[i + 1]
    return midsv_tags


def _complement_midsv(midsv_tags: list[str]) -> list[str]:
    comp = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N", "a": "t", "c": "g", "g": "c", "t": "a", "n": "n"}
    for i, cs in enumerate(midsv_tags):
        op = cs[0]
        if op == "*":
            midsv_tags[i] = op + comp[cs[1]] + comp[cs[2]]
        elif op == "+":
            midsv_tags[i] = "|".join(c[0] + comp[c[-1]] for c in cs.split("|"))
        else:  # Match or Deletion or N
            midsv_tags[i] = op + comp[cs[-1]]
        midsv_tags[i] = midsv_tags[i].replace("NN", "N")
        midsv_tags[i] = midsv_tags[i].replace("nn", "n")
    return midsv_tags


def revcomp_midsvs(midsv_tags: list[str]) -> list[str]:
    midsv_tags_reversed = _reverse_midsvs(midsv_tags)
    midsv_tags_realigned = _realign_insertion(midsv_tags_reversed)
    midsv_tags_revcomped = _complement_midsv(midsv_tags_realigned)
    return midsv_tags_revcomped


###########################################################
# convert midsv_tags to cstag
###########################################################


def _add_match_operator_to_n(midsv_tags: list[str]) -> list[str]:
    """Add "=" (match operator) to the sequences with "N"."""
    midsv_tags_add_match_op = []
    for cs in midsv_tags:
        if cs.startswith("+") and is_n_tag(cs.split("|")[-1]):
            cs_ins = cs.split("|")
            last_tag = cs_ins[-1]
            if not last_tag.startswith("="):
                last_tag = "=" + last_tag[-1]
            cs_ins[-1] = last_tag
            midsv_tags_add_match_op.append("|".join(cs_ins))
            continue

        if is_n_tag(cs):
            if cs.startswith("="):
                midsv_tags_add_match_op.append(cs)
            else:
                midsv_tags_add_match_op.append("=" + cs[-1])
            continue

        midsv_tags_add_match_op.append(cs)
    return midsv_tags_add_match_op


def _split_midsvs_by_delimiter(midsv_tags: list[str]) -> list[str]:
    midsv_tags_break = []
    for cs in midsv_tags:
        if cs.startswith("+"):
            midsv_tags_break.extend(cs.split("|"))
        else:
            midsv_tags_break.append(cs)
    return midsv_tags_break


def _combine_midsvs_by_prefix(midsv_tags: list[str]) -> list[str]:
    if not midsv_tags:
        return []

    combined_midsv_tags = []
    prev_prefix: str = midsv_tags[0][0]
    current_cs: list[str] = [midsv_tags[0][1:]]

    for item in midsv_tags[1:]:
        current_prefix, cs = item[0], item[1:]
        if prev_prefix == current_prefix:
            if current_prefix == "*":
                current_cs.append(item)
            else:
                current_cs.append(cs)
        else:
            combined_midsv_tags.append(prev_prefix + "".join(current_cs))
            prev_prefix = current_prefix
            current_cs = [cs]

    combined_midsv_tags.append(prev_prefix + "".join(current_cs))
    return combined_midsv_tags


def _standardize_case(midsv_tags: list[str]) -> list[str]:
    """Standardize the case of characters based on mutation types."""
    transformations = {"*": str.lower, "-": str.lower, "+": str.lower, "~": str.lower, "=": str.upper}

    transformed = []
    for cs in midsv_tags:
        prefix = cs[0]
        trans_func = transformations.get(prefix)
        transformed.append(prefix + trans_func(cs[1:]))

    return transformed


def convert_midsvs_to_cstag(midsv_tags: list[str]) -> str:
    midsv_tags_matched_op = _add_match_operator_to_n(midsv_tags)
    midsv_tags_splitted = _split_midsvs_by_delimiter(midsv_tags_matched_op)
    midsv_tags_combined = _combine_midsvs_by_prefix(midsv_tags_splitted)
    return "".join(_standardize_case(midsv_tags_combined))


###########################################################
# convert midsv_tags to sequence
###########################################################


def call_sequence(cons_percentage: list[dict[str, float]]) -> str:
    """convert position weight matrix (cons_pergentage) to sequence"""

    consensus_sequence = []
    for cons_per in cons_percentage:
        midsv_tags = max(cons_per, key=cons_per.get)
        cs_tag = convert_midsvs_to_cstag([midsv_tags])
        seq = cstag.to_sequence(cs_tag)
        consensus_sequence.append(seq)
    return "".join(consensus_sequence)
