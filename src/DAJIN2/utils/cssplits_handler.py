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


###########################################################
# reallocate_insertion_within_deletion
###########################################################

def highlight_sv_deletions(misdv_sv_allele: list[str]) -> list[str]:
    """Annotate start and end poisition at the SV deletions."""
    sv_deletions = []
    idx = 0
    while idx < len(misdv_sv_allele):
        midsv_tag = misdv_sv_allele[idx]

        if midsv_tag.startswith("-"):
            sv_deletions.append("!START_OF_DEL_ALLELE!")
            sv_deletions.append(misdv_sv_allele[idx])
            # Enclose consecutive deletions within a single span.
            while idx < len(misdv_sv_allele) - 1 and misdv_sv_allele[idx + 1].startswith("-"):
                sv_deletions.append(misdv_sv_allele[idx + 1])
                idx += 1
            sv_deletions.append("!END_OF_DEL_ALLELE!")
        # No SV
        else:
            sv_deletions.append(midsv_tag)

        idx += 1

    return sv_deletions


def embed_sv_deletions(midsv_consensus: list[str], midsv_sv_deletions: list[str]) -> list[str]:
    """
    Insert `!START_OF_DEL_ALLELE!` and `!END_OF_DEL_ALLELE!` into `midsv_consensus`.
    """
    idx = 0
    idx_sv_allele = 0

    midsv_consensus_with_deletion = []
    while idx_sv_allele < len(midsv_sv_deletions) and idx < len(midsv_consensus):
        tag_sv_allele = midsv_sv_deletions[idx_sv_allele]

        if tag_sv_allele.startswith("!"):
            midsv_consensus_with_deletion.append(tag_sv_allele)
            idx_sv_allele += 1
            continue

        if tag_sv_allele.startswith("-"):
            idx_sv_allele += 1
            continue

        tag_sample = midsv_consensus[idx]

        midsv_consensus_with_deletion.append(tag_sample)

        idx_sv_allele += 1
        idx += 1
    return midsv_consensus_with_deletion


def get_flanked_tags_by_deletions(midsv_tags: list[str]) -> tuple[list[list[str]], list[tuple[int, int]]]:
    """
    Extracts index and midsv tags flanked by large deletions.
    """
    flanked_tags = []
    flanked_indices = []
    i = 0
    while i < len(midsv_tags):
        tag = midsv_tags[i]
        flanked_tag = []

        # Flag indicating whether it is flanked until the next Del_Allele appears
        is_flanked = False

        if tag.endswith("!END_OF_DEL_ALLELE!") and i < len(midsv_tags) - 1:
            i += 1
            start_index = i
            while i < len(midsv_tags):
                if midsv_tags[i].startswith("!START_OF_DEL_ALLELE!"):
                    is_flanked = True
                    break
                flanked_tag.append(midsv_tags[i])
                i += 1

        if flanked_tag and is_flanked:
            flanked_tags.append(flanked_tag)
            end_index = i
            flanked_indices.append((start_index, end_index))

        i += 1

    return flanked_tags, flanked_indices


def has_consecutive_matches(tags: list[str], n: int = 10) -> bool:
    count = 0
    for tag in tags:
        if tag.startswith("="):
            count += 1
            if count >= n:
                return True
        else:
            count = 0
    return False


def get_last_index_of_del_allele(midsv_tags: list[str]) -> int:
    return next(i for i in range(len(midsv_tags) - 1, -1, -1) if midsv_tags[i].endswith("!END_OF_DEL_ALLELE!")) + 1


def convert_midsv_tags_to_insertion(midsv_tags: list[str]) -> str:
    """
    Converts MIDSV tags to insertion tags.
    e.g. ["=T", "=A", "+C|+C|=A", "-T", "=G", "=C"] -> "+T|+A|+C|+C|+A|+G|+C"
    """
    insertion_tags = []
    for tag in midsv_tags:
        if tag.startswith("-"):
            continue
        elif tag.startswith("+"):
            ins, last = "|".join(tag.split("|")[:-1]), tag.split("|")[-1]
            insertion_tags.append(ins)
            if last.startswith("-"):
                continue
            insertion_tags.append("+" + last[-1])
        else: # match or substitution
            insertion_tags.append("+" + tag[-1])
    return "|".join(insertion_tags)


def extract_deletion_tags(midsv_sv_deletions_highlighted: list[str]) -> list[list[str]]:
    deletion_tags = []
    for tag in midsv_sv_deletions_highlighted:
        if tag == "!START_OF_DEL_ALLELE!":
            deletion_tags.append([])
        if tag.startswith("-"):
            deletion_tags[-1].extend([tag])

    return deletion_tags


def reflect_deletions(midsv_consensus_with_deletion: list[str], deletion_tags: list[list[str]]) -> list[str]:
    """
    Reflect the Deletion when there is no Insertion.
    """
    for deletion_id, deletion_tag in enumerate(deletion_tags):
        idx_of_del = [i for i, tag in enumerate(midsv_consensus_with_deletion) if tag == "!END_OF_DEL_ALLELE!"][deletion_id]
        midsv_consensus_with_deletion[idx_of_del:idx_of_del] = deletion_tag

    return midsv_consensus_with_deletion


def handle_insertions_within_deletions(midsv_consensus_with_deletion: list[str], flanked_tags: list[list[str]]) -> list[str]:
    insertions = []
    i = 0

    while i < len(flanked_tags):
        flanked_tag = flanked_tags[i]

        # Get the index of the region enclosed by End and Start
        flank_start = [i for i, tag in enumerate(midsv_consensus_with_deletion) if tag == "!END_OF_DEL_ALLELE!"][i]
        flank_end = [i for i, tag in enumerate(midsv_consensus_with_deletion) if tag == "!START_OF_DEL_ALLELE!"][i+1]

        if has_consecutive_matches(flanked_tag) and insertions:
            # Reflect Insertions
            last_tag = midsv_consensus_with_deletion.pop(flank_start + 1)
            insertions = "|".join(["|".join(insertions), last_tag])
            midsv_consensus_with_deletion.insert(flank_start + 1, insertions)
            insertions = []

        elif not has_consecutive_matches(flanked_tag):
            # Collect Insertions
            insertions.append(convert_midsv_tags_to_insertion(flanked_tag))
            midsv_consensus_with_deletion[flank_start+1: flank_end] = ["!" for _ in range(flank_start+1, flank_end)]

        i += 1

    if insertions and "!END_OF_DEL_ALLELE!" in midsv_consensus_with_deletion:
        last_index = get_last_index_of_del_allele(midsv_consensus_with_deletion)
        last_tag = midsv_consensus_with_deletion.pop(last_index)
        insertions = "|".join(["|".join(insertions), last_tag])
        midsv_consensus_with_deletion.insert(last_index, insertions)

    # Remove ['!START_OF_DEL_ALLELE!', '!END_OF_DEL_ALLELE!'] tags
    return [m for m in midsv_consensus_with_deletion if not m.startswith("!")]

def reflect_sv_deletion_in_midsv(midsv_consensus: list[str], midsv_sv_deletions: list[str]) -> list[str]:
    """
    Since the mapping in minimap2 is local alignment, insertion bases within large deletions may be partially mapped to the reference genome and not detected as insertion bases. Therefore, update cssplits to detect insertions within large deletions as insertions.
    """
    midsv_sv_deletions_highlighted = highlight_sv_deletions(midsv_sv_deletions)
    midsv_consensus_with_deletion = embed_sv_deletions(midsv_consensus, midsv_sv_deletions_highlighted)

    # Append deletion tags to the midsv_consensus_with_deletion
    deletion_tags = extract_deletion_tags(midsv_sv_deletions_highlighted)
    midsv_consensus_with_deletion = reflect_deletions(midsv_consensus_with_deletion, deletion_tags)

    flanked_tags, _ = get_flanked_tags_by_deletions(midsv_consensus_with_deletion)

    return handle_insertions_within_deletions(midsv_consensus_with_deletion, flanked_tags)
