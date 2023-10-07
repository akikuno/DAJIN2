from __future__ import annotations


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
# convert cssplits to cstag
###########################################################


def add_match_operator_to_n(cssplits: list[str]) -> list[str]:
    """Add "=" (match operator) to the sequences that start with "N"."""
    return ["=" + seq if seq.startswith("N") else seq for seq in cssplits]


def concatenate_cssplits(cssplits: list[str]) -> str:
    """Concatenate list of sequences based on certain variants."""
    concatenated = [cssplits[0]]

    for prev, current in zip(cssplits, cssplits[1:]):
        if prev[0] == current[0] and current[0] in {"=", "-"}:
            concatenated.append(current[1])
        elif current[0] == "+":
            cs_ins = current.replace("+", "").split("|")
            concatenated.append("+" + "".join(cs_ins[:-1]) + cs_ins[-1])
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
