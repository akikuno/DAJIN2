from __future__ import annotations

import re

###############################################################################
# HTML Builder
###############################################################################

HTML_HEADER = """<!DOCTYPE html>
<html lang="en">
<meta charset="utf-8">

<head>
    <style>
        body {
            width: 95%;
            margin: 0 auto;
            box-sizing: border-box;
            font-family: Consolas, 'Courier New', monospace;
            color: #333;
        }

        h1 {
            padding: 0.1em 0;
            border-top: solid 3px #333;
            border-bottom: solid 3px #333;
        }

        .p_seq {
            color: #585858;
            word-wrap: break-word;
            letter-spacing: 0.15em;
        }

        .p_legend {
            word-wrap: break-word;
        }

        .p_legend span {
            margin: 5px;
        }

        .p_legend br {
            line-height: 2em;
        }

        .Ins {
            font-weight: bold;
            border: 0.1em solid;
            background-color: #ee827c;
        }

        .Del {
            font-weight: bold;
            border: 0.1em solid;
            background-color: #a0d8ef;
        }

        .Sub {
            font-weight: bold;
            border: 0.1em solid;
            background-color: #98d98e;
        }

        .Inv {
            font-weight: bold;
            border: 0.1em solid;
            background-color: #CAB5FF;
        }

        .Splice {
            font-weight: bold;
            border: 0.1em solid;
            background-color: #f8e58c;
        }

        .Unknown {
            font-weight: bold;
            border: 0.1em solid;
            background-color: #c0c6c9;
        }

        .Ins_Allele {
            font-weight: bold;
            text-decoration: underline;
            text-decoration-color: #ED5F5F;
            text-decoration-thickness: 3px;
        }

        .Del_Allele {
            font-weight: bold;
            text-decoration: underline;
            text-decoration-color: #48C0F0;
            text-decoration-thickness: 3px;
        }

        .Inv_Allele {
            font-weight: bold;
            text-decoration: underline;
            text-decoration-color: #7F4DFF;
            text-decoration-thickness: 3px;
        }
    </style>
</head>

<body>
"""

HTML_LEGEND = """
<p class="p_legend">
    <span class="Ins">Insertion</span>
    <span class="Del">Deletion</span>
    <span class="Sub">Substitution</span>
    <span class="Splice">Splice</span>
    <span class="Inv">Inversion</span>
    <span class="Unknown">Unknown</span>
    <br>
    <span class="Ins_Allele">Insertion allele</span>
    <span class="Del_Allele">Deletion allele</span>
    <span class="Inv_Allele">Inversion allele</span>
</p>
<hr>
"""

HTML_FOOTER = """
</body>
</html>
"""


###############################################################################
# Highlight SV regions
###############################################################################

def highlight_sv_regions(misdv_sv_allele: list[str]) -> list[str]:
    """Add an HTML class for a colored underline to regions corresponding to SVs (indels, inversions)."""
    html_sv_allele = []
    idx = 0
    while idx < len(misdv_sv_allele):
        midsv_tag = misdv_sv_allele[idx]
        
        # Insertion
        if midsv_tag.startswith("+"):
            html_sv_allele.append("<span class='Ins_Allele'>")
            insertions, last_tag = midsv_tag.split("|")[:-1], midsv_tag.split("|")[-1]
            ins_seq = ["=" + ins[1] for ins in insertions]
            html_sv_allele.extend(ins_seq)
            html_sv_allele.append("</span>")
            # If a deletion occurs immediately after an insertion, close the insertion span and open the deletion span.
            if last_tag.startswith("-"):
                html_sv_allele.append("<span class='Del_Allele'>")
                html_sv_allele.append(last_tag)
                html_sv_allele.append("</span><!-- END_OF_DEL_ALLELE -->")
            else:  # match or substitution
                html_sv_allele.append(last_tag)
        
        # Deletion
        elif midsv_tag.startswith("-"):
            html_sv_allele.append("<span class='Del_Allele'>")
            html_sv_allele.append(misdv_sv_allele[idx])
            # Enclose consecutive deletions within a single span.
            while idx < len(misdv_sv_allele) - 1 and misdv_sv_allele[idx + 1].startswith("-"):
                html_sv_allele.append(misdv_sv_allele[idx + 1])
                idx += 1
            html_sv_allele.append("</span><!-- END_OF_DEL_ALLELE -->")

        # Inversion
        elif midsv_tag.islower():
            html_sv_allele.append("<span class='Inv_Allele'>")
            # Enclose consecutive inversion within a single span.
            html_sv_allele.append(midsv_tag)
            while idx < len(misdv_sv_allele) - 1 and misdv_sv_allele[idx + 1].islower():
                html_sv_allele.append(misdv_sv_allele[idx + 1])
                idx += 1
            html_sv_allele.append("</span>")

        # No SV
        else:
            html_sv_allele.append(midsv_tag)

        idx += 1

    return html_sv_allele

###############################################################################
# Embed mutations to HTML
###############################################################################

def embed_mutations_to_html(html_sv_allele: list[str], midsv_consensus: list[str]) -> list[str]:
    idx = 0
    idx_sv_allele = 0
    html_parts = ["<p class='p_seq'>"]

    while idx_sv_allele < len(html_sv_allele):
        tag_sv_allele = html_sv_allele[idx_sv_allele]

        if idx >= len(midsv_consensus):
            break

        if tag_sv_allele.startswith("<"):
            html_parts.append(tag_sv_allele)
            idx_sv_allele += 1
            continue

        if tag_sv_allele.startswith("-"):
            html_parts.append(tag_sv_allele[-1])
            idx_sv_allele += 1
            continue

        tag_sample = midsv_consensus[idx]

        # Insertion
        if tag_sample.startswith("+"):
            insertions, last_tag = tag_sample.split("|")[:-1], tag_sample.split("|")[-1]
            insertions = "".join(ins[-1] for ins in insertions)
            html_parts.append("<span class='Ins'>" + "".join(insertions) + "</span>")

            if last_tag.startswith("-"):
                html_parts.append("<span class='Del'>" + last_tag[-1] + "</span>")
            else:  # match or substitution
                html_parts.append(last_tag[-1])

        # Substitution
        elif tag_sample.startswith("*"):
            substitutions = [tag_sample[-1]]
            while idx < len(midsv_consensus) - 1 and midsv_consensus[idx + 1].startswith("*"):
                substitutions.append(midsv_consensus[idx + 1][-1])
                idx_sv_allele += 1
                idx += 1
            html_parts.append("<span class='Sub'>" + "".join(substitutions) + "</span>")

        # Deletion
        elif tag_sample.startswith("-"):
            deletions = [tag_sample[-1]]
            while idx < len(midsv_consensus) - 1 and midsv_consensus[idx + 1].startswith("-"):
                deletions.append(midsv_consensus[idx + 1][-1])
                idx_sv_allele += 1
                idx += 1
            html_parts.append("<span class='Del'>" + "".join(deletions) + "</span>")

        # Othres
        else:
            html_parts.append(tag_sample[-1])

        idx_sv_allele += 1
        idx += 1

    html_parts.append("</p>")

    return html_parts


###############################################################################
# Handle insertions
###############################################################################

def get_flanked_tags_by_deletions(html_parts: list[str]) -> tuple[list[list[str]], list[tuple[int, int]]]:
    """
    Extracts regions flanked by large deletions from an HTML parts list.

    This function scans through a list of HTML parts and identifies segments that 
    are enclosed between the comment `<!-- END_OF_DEL_ALLELE -->` and the next occurrence 
    of `<span class='Del_Allele'>`. These regions are extracted along with their start 
    and end indices.
    """
    flanked_tags = []
    flanked_indices = []
    i = 0
    while i < len(html_parts):
        tag = html_parts[i]
        flanked_tag = []

        # Flag indicating whether it is flanked until the next Del_Allele appears
        is_flanked = False

        if tag.endswith("<!-- END_OF_DEL_ALLELE -->"):
            i += 1
            start_index = i
            while i < len(html_parts):
                if html_parts[i].startswith("<span class='Del_Allele'>"):
                    is_flanked = True
                    break
                flanked_tag.append(html_parts[i])
                i += 1

        if flanked_tag and is_flanked:
            flanked_tags.append(flanked_tag)
            end_index = i
            flanked_indices.append((start_index, end_index))

        i += 1

    return flanked_tags, flanked_indices


def get_sequence_from_flanked_tags(flanked_tag: list[str]) -> str:
    flanked_tag_remove_del = [tag for tag in flanked_tag if not tag.startswith("<span class='Del'>")]
    sequence = re.findall(r'[ACTGN]', "".join(flanked_tag_remove_del))
    return sequence


def has_consecutive_matches(flanked_tag: list[str], n: int = 10) -> bool:
    count_match = 0
    for tag in flanked_tag:
        if tag.startswith("<"):
            count_match = 0
        else:
            count_match += 1
        
    if count_match >= n:
        return True
    
    return False


def get_last_index_of_del_allele(html_parts: list[str]) -> int:
    return next(i for i in range(len(html_parts) - 1, -1, -1) if html_parts[i].endswith("<!-- END_OF_DEL_ALLELE -->")) + 1



def handle_insertions(html_parts: list[str], flanked_tags: list[list[str]], flanked_indices: list[tuple[int, int]], sequences_within_deletion: list[str]) -> list[str]:
    insertions = []
    previous_end = None
    i = 0
    index_shift = 0
    while i < len(flanked_tags):
        flanked_tag = flanked_tags[i]
        start, end = flanked_indices[i]

        start += index_shift
        end += index_shift

        sequence_within_deletion =  sequences_within_deletion[i]

        if has_consecutive_matches(flanked_tag) and insertions:
            insertions = "".join(insertions)
            html_parts.insert(previous_end + 2, f"<span class='Ins'>{insertions}</span>")
            index_shift += 1
            insertions = []

        elif not has_consecutive_matches(flanked_tag):
            insertions.append("".join(get_sequence_from_flanked_tags(flanked_tag)))

            # print(html_parts[start-1], html_parts[end])

            html_parts[start-1: end+1] = sequence_within_deletion
            index_shift -= 2

        # print(html_parts)
        previous_end = end
        i += 1

    if insertions:
        last_index = get_last_index_of_del_allele(html_parts)
        insertions = "".join(insertions)
        html_parts.insert(last_index, f"<span class='Ins'>{insertions}</span>")

    return html_parts

###############################################################################
# main
###############################################################################


def to_html(midsv_sv_allele: list[str], midsv_consensus: list[str:], description: str = "") -> str:
    """Output HTML string showing a sequence with mutations colored"""
    description_str = f"<h1>{description}</h1>" if description else ""
    html_sv_allele = highlight_sv_regions(midsv_sv_allele)
    html_parts = embed_mutations_to_html(html_sv_allele, midsv_consensus)

    # Handle insertions within large deletions
    flanked_tags_sv_allele, _ = get_flanked_tags_by_deletions(html_sv_allele)
    sequences_within_deletion = [get_sequence_from_flanked_tags(tag) for tag in flanked_tags_sv_allele]
    flanked_tags, flanked_indices = get_flanked_tags_by_deletions(html_parts)
    html_parts = "".join(handle_insertions(html_parts, flanked_tags, flanked_indices, sequences_within_deletion))

    return "\n".join(
        [
            HTML_HEADER,
            description_str,
            HTML_LEGEND,
            html_parts,
            HTML_FOOTER,
        ]
    )
