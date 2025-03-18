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

        .Unknown {
            font-weight: bold;
            border: 0.1em solid;
            background-color: #c0c6c9;
        }

        .Inv {
            font-weight: bold;
            text-decoration: underline;
            text-decoration-color: #7F4DFF;
            text-decoration-thickness: 3px;
            text-underline-offset: auto;
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
    <span class="Unknown">Unknown</span>
    <span class="Inv">Inversion</span>
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


def get_inversion_range(midsv_sv_allele) -> list[tuple[int, int]]:
    inversion_range = []
    prev_inversion = False
    inversion_start = -1
    for i, tag_sv in enumerate(midsv_sv_allele):
        if tag_sv.islower() and prev_inversion is False:
            inversion_start = i
            prev_inversion = True
        elif not tag_sv.islower() and prev_inversion is True:
            inversion_range.append((inversion_start, i - 1))
            inversion_start = -1
            prev_inversion = False

        # Last tag is inversion
        if i == len(midsv_sv_allele) - 1 and prev_inversion is True:
            inversion_range.append((inversion_start, i))

    return inversion_range


def get_insertion_range(midsv_sv_allele) -> list[tuple[int, int]]:
    insertion_range = []
    idx = 0
    for i, tag_sv in enumerate(midsv_sv_allele):
        if tag_sv.startswith("+"):
            insertion_start = i + idx
            insertion_end = insertion_start + tag_sv.count("|") - 1
            insertion_range.append((insertion_start, insertion_end))
            idx += insertion_end
    return insertion_range


def add_sv_class(cons_midsv_tag, sv_range, class_name: str) -> list[str]:
    sv_range = iter(sv_range)
    highlight_sv_allele = []

    sv_start = -1
    sv_end = -1
    for i, tag_sample in enumerate(cons_midsv_tag):
        if sv_start == -1 or i == sv_end + 1:
            try:
                sv_start, sv_end = next(sv_range)
            except StopIteration:
                pass
        if i == sv_start == sv_end:
            highlight_sv_allele.append(f"<span class='{class_name}'>{tag_sample}</span>")
        elif i == sv_start:
            highlight_sv_allele.append(f"<span class='{class_name}'>{tag_sample}")
        elif i == sv_end:
            highlight_sv_allele.append(f"{tag_sample}</span>")
        else:
            highlight_sv_allele.append(tag_sample)

    return highlight_sv_allele


def append_inversion_allele(midsv_sv_allele, cons_midsv_tag) -> list[str]:
    inversion_range = get_inversion_range(midsv_sv_allele)
    return add_sv_class(cons_midsv_tag, inversion_range, "Inv")


def append_insertion_allele(midsv_sv_allele, cons_midsv_tag) -> list[str]:
    insertion_range = get_insertion_range(midsv_sv_allele)
    return add_sv_class(cons_midsv_tag, insertion_range, "Ins_Allele")


def split_html_tags(highlight_sv_allele: list[str]) -> list[str]:
    """
    Splits a list of strings by HTML tags.

    Example:
        ["<span>=A</span>", "<span>-C", "=C", "</span>"]
        -> ["<span>", "=A", "</span>", "<span>", "-C", "=C", "</span>"]
    """
    pattern = r"(<[^>]+>)"  # Regular expression to capture HTML tags
    return [part for text in highlight_sv_allele for part in re.split(pattern, text) if part]


###############################################################################
# Embed mutations to HTML
###############################################################################


def embed_mutations_to_html(highlight_sv_allele: list[str]) -> list[str]:
    idx = 0
    html_parts = ["<p class='p_seq'>"]

    while idx < len(highlight_sv_allele):
        tag_sample = highlight_sv_allele[idx]

        if tag_sample.startswith("<"):  # HTML tag
            html_parts.append(tag_sample)
            idx += 1
            continue

        # Insertion
        if tag_sample.startswith("+"):
            insertions, last_tag = tag_sample.split("|")[:-1], tag_sample.split("|")[-1]
            insertions = "".join(ins[-1] for ins in insertions)
            html_parts.append("<span class='Ins'>" + "".join(insertions) + "</span>" + last_tag[-1])
            idx += 1
            continue

        # Substitution
        if tag_sample.startswith("*"):
            substitutions = []
            while idx < len(highlight_sv_allele) and highlight_sv_allele[idx].startswith("*"):
                substitutions.append(highlight_sv_allele[idx][-1])
                idx += 1
            html_parts.append("<span class='Sub'>" + "".join(substitutions) + "</span>")
            continue

        # Deletion
        if tag_sample.startswith("-"):
            deletions = []
            while idx < len(highlight_sv_allele) and highlight_sv_allele[idx].startswith("-"):
                deletions.append(highlight_sv_allele[idx][-1])
                idx += 1
            html_parts.append("<span class='Del'>" + "".join(deletions) + "</span>")
            continue

        # Inversion: Just convert to upper case
        if tag_sample.islower() and not tag_sample.startswith("<"):  # avoid HTML span tag
            inversions = []
            while (
                idx < len(highlight_sv_allele)
                and highlight_sv_allele[idx].islower()
                and not highlight_sv_allele[idx].startswith("<")
            ):
                tag_inversion = highlight_sv_allele[idx].upper()
                if tag_inversion.startswith("*"):  # Substitution
                    inversions.append("<span class='Sub'>" + tag_inversion[-1] + "</span>")
                elif tag_inversion.startswith("-"):  # Deletion
                    inversions.append("<span class='Del'>" + tag_inversion[-1] + "</span>")
                elif tag_inversion.startswith("+"):  # Insertion
                    tag_insertion, tag_last = tag_inversion.split("|")[:-1], tag_inversion.split("|")[-1]
                    insertion = "".join(ins[-1] for ins in tag_insertion)
                    inversions.append("<span class='Ins'>" + "".join(insertion) + "</span>" + tag_last[-1])
                else:
                    inversions.append(tag_inversion[-1])
                idx += 1
            html_parts.append("".join(inversions))
            continue

        # Othres
        html_parts.append(tag_sample[-1])
        idx += 1

    html_parts.append("</p>")

    return html_parts


###############################################################################
# main
###############################################################################


def to_html(
    midsv_sv_allele: list[str], cons_midsv_tag: list[str:], allele: str, is_sv_allele: bool, description: str = ""
) -> str:
    """Output HTML string showing a sequence with mutations colored"""
    description_str = f"<h1>{description}</h1>" if description else ""
    if allele.startswith("insertion") and is_sv_allele:
        highlight_sv_allele = append_insertion_allele(midsv_sv_allele, cons_midsv_tag)
    elif allele.startswith("inversion") and is_sv_allele:
        highlight_sv_allele = append_inversion_allele(midsv_sv_allele, cons_midsv_tag)
    else:
        highlight_sv_allele = cons_midsv_tag

    highlight_sv_allele = split_html_tags(highlight_sv_allele)

    html_parts = embed_mutations_to_html(highlight_sv_allele)
    html_parts_str = "".join(html_parts)

    return "\n".join(
        [
            HTML_HEADER,
            description_str,
            HTML_LEGEND,
            html_parts_str,
            HTML_FOOTER,
        ]
    )
