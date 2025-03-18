from __future__ import annotations

###############################################################################
# annotate_sv_allele
###############################################################################


def annotate_sv_allele(cons_midsv_tag: list[str], sv_midsv_tag: list[str]) -> list[str]:
    """Annotate the SV allele to consenses midsv tag."""
    cons_midsv_tag_with_sv_tag = []
    idx = 0
    num_deletion = 0
    num_insertion = 0
    while idx < len(sv_midsv_tag) and idx - num_deletion + num_insertion < len(cons_midsv_tag):
        sv_tag = sv_midsv_tag[idx]

        # Deletion
        if sv_tag.startswith("-"):
            cons_midsv_tag_with_sv_tag.append(sv_tag)
            # Enclose consecutive alleles within a single span.
            num_deletion += 1
            idx += 1
            while idx < len(sv_midsv_tag) and sv_midsv_tag[idx].startswith("-"):
                cons_midsv_tag_with_sv_tag.append(sv_midsv_tag[idx])
                num_deletion += 1
                idx += 1

        # Inversion
        elif sv_tag.islower():
            cons_midsv_tag_with_sv_tag.append(cons_midsv_tag[idx].lower())
            idx += 1

        # Insertion
        elif sv_tag.startswith("+"):
            idx_insertion = idx + num_insertion
            end_insertion = idx_insertion + sv_tag.count("|")
            insertions = []
            while idx_insertion < len(cons_midsv_tag) and idx_insertion < end_insertion:
                cons_tag = cons_midsv_tag[idx_insertion]

                if cons_tag.startswith("=") or cons_tag.startswith("*"):
                    insertions.append(f"+{cons_tag[-1]}")

                elif cons_tag.startswith("+"):
                    ins_tag = "|".join(cons_tag.split("|")[:-1])
                    last_tag = cons_tag.split("|")[-1]
                    if last_tag.startswith("-"):
                        last_tag = ""
                    else:
                        last_tag = f"+{last_tag[-1]}"
                    insertions.append("|".join([ins_tag, last_tag]))
                elif cons_tag.startswith("-"):
                    pass

                idx_insertion += 1

            # last_tagが欠失の場合にスキップする
            while cons_midsv_tag[idx_insertion].startswith("-"):
                num_insertion += 1
                idx_insertion += 1

            insertions = f"{'|'.join(insertions)}|{cons_midsv_tag[idx_insertion]}"
            cons_midsv_tag_with_sv_tag.append(insertions)

            num_insertion += sv_tag.count("|")
            idx += 1
        # No SV
        else:
            cons_midsv_tag_with_sv_tag.append(cons_midsv_tag[idx - num_deletion + num_insertion])
            idx += 1
    return cons_midsv_tag_with_sv_tag


###############################################################################
# annotate_insertions_within_deletion
""" Strategy
- Extract regions flanked by deletions.
- Determine whether the regions are insertions.
- Convert tags to insertion format.
- Merge consecutive insertions.
- Annotate insertions in the consensus midsv tag list.
"""
###############################################################################


def extract_flanked_ranges(sv_midsv_tag: list[str]) -> list[tuple[int, int]]:
    """Extract start and end indices of regions flanked by deletions ("-")."""
    prev_deletion = False
    start = -1
    end = -1
    flanked_ranges = []
    i = 0

    while i < len(sv_midsv_tag):
        sv_tag = sv_midsv_tag[i]
        cur_deletion = sv_tag.startswith("-")

        if prev_deletion and not cur_deletion:
            start = i
        elif not prev_deletion and cur_deletion and start != -1:
            end = i - 1

        if start != -1 and end != -1:
            flanked_ranges.append((start, end))
            start = -1
            end = -1

        prev_deletion = cur_deletion
        i += 1

    # Append sentinel empty value to match the number of deletions, as location candidates are one more than tag_flanked.
    flanked_ranges.append(())

    return flanked_ranges


def get_end_of_deletion_indices(sv_midsv_tag: list[str]) -> list[int]:
    """Retrieve the indices marking the end of deletion regions."""
    index_end_of_deletion = []
    is_prev_del = False

    for i, tag in enumerate(sv_midsv_tag):
        if tag.startswith("-"):
            is_prev_del = True
        elif is_prev_del:
            index_end_of_deletion.append(i)
            is_prev_del = False

    return index_end_of_deletion


def extract_tags_of_flanked_ranges(cons_midsv_tag: list[str], index_flanked: list[tuple[int, int]]) -> list[str]:
    """Extract tags of regions flanked by deletions ("-")."""
    tag_flanked = []
    for indices in index_flanked:
        if indices:
            tag_flanked.append([cons_midsv_tag[i] for i in range(indices[0], indices[1] + 1)])
    return tag_flanked


def detect_insertions(tag_flanked: list[list[str]], n_match: int = 10) -> list[bool]:
    """Determine whether the regions flanked by deletions are insertion sequences."""
    count_consective_match = 0
    is_insertions = []
    for tags in tag_flanked:
        is_insertion = True
        for tag in tags:
            if tag.startswith("="):
                count_consective_match += 1
            else:
                count_consective_match = 0
            if count_consective_match >= n_match:
                is_insertion = False
                break
        is_insertions.append(is_insertion)

    # Append sentinel value to match the number of deletions, as location candidates are one more than tag_flanked.
    is_insertions.append(False)

    return is_insertions


def convert_tags_to_insertion(tag_flanked: list[list[str]], is_insertions: list[bool]) -> list[list[str]]:
    """Convert tags to insertion format based on the given insertion flags."""
    tag_insertion = []

    for tags, is_insertion in zip(tag_flanked, is_insertions):
        if is_insertion:
            new_tags = [f"+{tag[-1]}" for tag in tags if not tag.startswith("-")]
        else:
            new_tags = tags  # Keep original tags unchanged
        tag_insertion.append(new_tags)

    return tag_insertion


def merge_inserted_sequences(tag_insertion: list[list[str]], is_insertions: list[bool]) -> list[list[str] | None]:
    """Merge insertion regions flanked by deletions."""
    tag_insertion_merged = []
    tmp_insertions = []

    i = 0
    while i < len(tag_insertion):
        if is_insertions[i]:
            tag_insertion_merged.append([])  # Mark the start of an insertion
            tmp_insertions += tag_insertion[i]

            # Merge consecutive insertions
            if i + 1 < len(is_insertions) and is_insertions[i + 1]:
                tmp_insertions += tag_insertion[i + 1]
                tag_insertion[i + 1] = []  # Mark as processed
            else:
                tag_insertion_merged.append(tmp_insertions)
                tmp_insertions = []

            tag_insertion[i] = []  # Mark current as processed
        else:
            tag_insertion_merged.append(None)  # Non-insertion region remains unchanged

        i += 1

    return tag_insertion_merged


def annotate_insertion(
    cons_midsv_tag: list[str],
    sv_midsv_tag: list[str],
    tag_insertion_merged: list[list[str] | None],
    flanked_ranges: list[list[int]],
    index_end_of_deletion: list[int],
) -> list[str]:
    """
    Annotate insertions in the consensus midsv tag list.
    """
    for tag_insertion, idx_flanked, idx_del in zip(tag_insertion_merged, flanked_ranges, index_end_of_deletion):
        if tag_insertion is None:
            continue

        if not tag_insertion and idx_flanked:
            # Replace the flanked region with deletion tags
            sv_deletions = [tag.replace("=", "-") for tag in sv_midsv_tag[idx_flanked[0] : idx_flanked[1] + 1]]
            cons_midsv_tag[idx_flanked[0] : idx_flanked[1] + 1] = sv_deletions
        else:
            i = idx_del
            # Find the first non-deletion tag
            while i < len(cons_midsv_tag) and cons_midsv_tag[i].startswith("-"):
                i += 1

            # Append the existing tag and join with "|"
            if i < len(cons_midsv_tag):  # Ensure within bounds
                tag_insertion.append(cons_midsv_tag[i])
                cons_midsv_tag[i] = "|".join(tag_insertion)

    return cons_midsv_tag


###############################################################################
# main
###############################################################################


def annotate_insertions_within_deletion(cons_midsv_tag: list[str], sv_midsv_tag: list[str]) -> list[str]:
    """Reflect SV deletions in the consensus midsv tag list."""
    flanked_ranges = extract_flanked_ranges(sv_midsv_tag)
    index_end_of_deletion = get_end_of_deletion_indices(sv_midsv_tag)
    tag_flanked = extract_tags_of_flanked_ranges(cons_midsv_tag, flanked_ranges)
    is_insertions = detect_insertions(tag_flanked)
    tag_insertion = convert_tags_to_insertion(tag_flanked, is_insertions)
    tag_insertion_merged = merge_inserted_sequences(tag_insertion, is_insertions)
    return annotate_insertion(
        cons_midsv_tag, sv_midsv_tag, tag_insertion_merged, flanked_ranges, index_end_of_deletion
    )
