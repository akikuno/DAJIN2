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
    <span class="Inv">Inversion</span>
    <span class="Unknown">Unknown</span>
    <br>
    <span>SV:</span>
    <span class="Ins_Allele">Insertion allele</span>
    <span class="Del_Allele">Insertion allele</span>
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
            html_sv_allele.append(last_tag)

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

    while idx_sv_allele < len(html_sv_allele) and idx < len(midsv_consensus):

        tag_sv_allele = html_sv_allele[idx_sv_allele]
        if tag_sv_allele.startswith("<"):
            html_parts.append(tag_sv_allele)
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
# main
###############################################################################


def to_html(midsv_sv_allele: list[str], midsv_consensus: list[str:], description: str = "") -> str:
    """Output HTML string showing a sequence with mutations colored"""
    description_str = f"<h1>{description}</h1>" if description else ""
    html_sv_allele = highlight_sv_regions(midsv_sv_allele)
    html_parts = embed_mutations_to_html(html_sv_allele, midsv_consensus)
    html_parts = "".join(html_parts)

    return "\n".join(
        [
            HTML_HEADER,
            description_str,
            HTML_LEGEND,
            html_parts,
            HTML_FOOTER,
        ]
    )
