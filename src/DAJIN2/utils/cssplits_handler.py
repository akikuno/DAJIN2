from __future__ import annotations

import re

import cstag
import numpy as np
import ruptures as rpt


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


def call_sequence(cons_percentage: list[dict[str, float]], sep: str = "") -> str:
    """convert position weight matrix (cons_pergentage) to sequence"""

    consensus_sequence = []
    for cons_per in cons_percentage:
        cssplits = max(cons_per, key=cons_per.get)
        cs_tag = convert_cssplits_to_cstag([cssplits])
        seq = cstag.to_sequence(cs_tag)
        consensus_sequence.append(seq)
    return f"{sep}".join(consensus_sequence)


###########################################################
# reallocate_insertion_within_deletion
###########################################################


def _extract_candidate_index_of_large_deletions(
    cssplits: list[str], bin_size: int = 500, percentage: int = 50
) -> list[dict[str, int]]:
    deletion_counts = []
    for start in range(0, len(cssplits), bin_size // 10):
        end = start + bin_size
        counts_deletion = "".join(cssplits[start:end]).count("-")
        deletion_counts.append((start, counts_deletion))

    range_of_large_deletions = []
    on_the_peak = False
    start = -1
    end = -1
    for i, count in deletion_counts:
        del_percentage = count / bin_size * 100
        if del_percentage > percentage and on_the_peak is False:
            start = i
            on_the_peak = True
        if del_percentage <= percentage and on_the_peak is True:
            end = i
            on_the_peak = False
        if start != -1 and end != -1:
            range_of_large_deletions.append({"start": start, "end": end})
            start = -1
            end = -1

    return range_of_large_deletions


def _extract_break_points_of_large_deletions(
    cssplits: list[str], range_of_large_deletions: list[dict[str, int]], bin_size: int = 500
):
    cssplits_match_bool = np.array([1 if cs.startswith("=") else 0 for cs in cssplits], dtype=bool)

    break_points = []
    for range_of_deletion in range_of_large_deletions:
        break_point = {}
        for key, index in range_of_deletion.items():
            window_start = index
            window_end = index + bin_size

            # Find breakpoints(bps). First, predict approximate breakpoints using the signal of the window.
            signal = cssplits_match_bool[window_start:window_end]
            bps = rpt.Window(model="ar").fit(signal).predict(n_bkps=10)[:-1]

            # Aling breakpoints to the original index
            if bps:
                approximate_bp: int = bps[0] if key == "start" else bps[-1]
            else:
                approximate_bp: int = 0 if key == "start" else bin_size

            # Find the exact breakpoint because the approximate breakpoint may not be accurate
            signal_subset = signal[max(0, approximate_bp - 50) : min(approximate_bp + 50, len(signal))]
            bp = -1
            for i, sig in enumerate(signal_subset):
                if key == "start" and not sig:
                    bp = i + max(0, approximate_bp - 50) + window_start
                    break
                if key == "end" and sig:
                    bp = i + max(0, approximate_bp - 50) + window_start - 1
                    break
            if bp == -1:
                bp = approximate_bp

            # Store the breakpoint
            break_point[key] = bp

        if break_point["start"] == break_point["end"]:
            continue

        break_points.append(break_point)

    return break_points


def _convert_break_points_to_index(break_points: list[dict[str, int]]) -> list[int]:
    index_of_large_deletions = []
    for break_point in break_points:
        start = break_point["start"]
        end = break_point["end"]
        index_of_large_deletions += list(range(start, end + 1))

    return index_of_large_deletions


def _find_matched_indexes(cssplits: list[str], index_of_large_deletions: list[int]) -> list[int]:
    matched_index = []
    count_matches = 0
    start_match = -1

    index_of_large_deletions.sort()
    for i in index_of_large_deletions:
        if cssplits[i].startswith("="):
            if start_match == -1:
                start_match = i
            count_matches += 1
        else:
            if count_matches >= 10:
                matched_index += list(range(start_match, i))
            count_matches = 0
            start_match = -1

    return matched_index


def _remove_matched_indexes(index_of_large_deletions: list[int], matched_index: list[int]) -> set[int]:
    return set(index_of_large_deletions) - set(matched_index)


def _get_index_of_large_deletions(cssplits: list[str], bin_size: int = 500, percentage: int = 50) -> set[int]:
    range_of_large_deletions = _extract_candidate_index_of_large_deletions(cssplits, bin_size, percentage)
    break_points = _extract_break_points_of_large_deletions(cssplits, range_of_large_deletions, bin_size)

    index_of_large_deletions = _convert_break_points_to_index(break_points)
    matched_index = _find_matched_indexes(cssplits, index_of_large_deletions)
    return _remove_matched_indexes(index_of_large_deletions, matched_index)


def _adjust_cs_insertion(cs: str) -> str:
    if cs.startswith("-"):
        return None
    elif cs.startswith("+"):
        ins = "|".join(cs.split("|")[:-1])
        end_cs = cs.split("|")[-1][-1]
        cs_insertion = f"{ins}|+{end_cs}|"
    else:
        cs_insertion = f"+{cs[-1]}|"
    return cs_insertion


def reallocate_insertion_within_deletion(cssplits: list[str], bin_size: int = 500, percentage: int = 50) -> list[str]:
    """
    Since the mapping in minimap2 is local alignment, insertion bases within large deletions may be partially mapped to the reference genome and not detected as insertion bases. Therefore, update cssplits to detect insertions within large deletions as insertions.
    """
    cssplits_updated = cssplits.copy()

    index_of_large_deletions: set[int] = _get_index_of_large_deletions(
        cssplits, bin_size=bin_size, percentage=percentage
    )

    insertion_within_deletion = []
    for i, cs in enumerate(cssplits):
        if i in index_of_large_deletions:
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

            insertion_within_deletion.append(_adjust_cs_insertion(cs))
        else:
            if not insertion_within_deletion:
                continue
            cssplits_updated[i] = "".join(insertion_within_deletion) + cs
            insertion_within_deletion = []

    return cssplits_updated
