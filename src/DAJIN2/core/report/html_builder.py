HTML_HEADER = """<!DOCTYPE html>
    <html>
    <head>
    <style>
    h1 {
    font-family: Consolas, 'Courier New', monospace;
    color: #333;
    padding: 0.1em 0;
    border-top: solid 3px #333;
    border-bottom: solid 3px #333;
    }
    .p_seq {
    font-family: Consolas, 'Courier New', monospace;
    color: #585858;
    word-wrap: break-word;
    letter-spacing: 0.15em;
    }
    .p_legend {
    font-family: Consolas, 'Courier New', monospace;
    color: #585858;
    word-wrap: break-word;
    }

    <!-- Annotating mutations -->
    .Ins {
    color: #333;
    font-weight: bold;
    border: 0.1em solid;
    background-color: #ee827c;
    }
    .Del {
    color: #333;
    font-weight: bold;
    border: 0.1em solid;
    background-color: #a0d8ef;
    }
    .Sub {
    color: #333;
    font-weight: bold;
    border: 0.1em solid;
    background-color: #98d98e;
    }
    .Inv {
    color: #333;
    font-weight: bold;
    border: 0.1em solid;
    background-color: #CAB5FF;
    }
    .Splice {
    color: #333;
    font-weight: bold;
    border: 0.1em solid;
    background-color: #f8e58c;
    }
    .Unknown {
    color: #333;
    font-weight: bold;
    border: 0.1em solid;
    background-color: #c0c6c9;
    }
    
    <!-- Annotating SV to the consensus allele -->
    .Ins_Allele {
    color: #333;
    font-weight: bold;
    text-decoration:underline;
    text-decoration-color: #ED5F5F;
    }
    .Del_Allele {
    color: #333;
    font-weight: bold;
    text-decoration:underline;
    text-decoration-color: #48C0F0;
    }
    .Inv_Allele {
    color: #333;
    font-weight: bold;
    text-decoration:underline;
    text-decoration-color: #7F4DFF;
    }
    </style>
    </head>
    <body>
"""

HTML_LEGEND = """
<p class = "p_legend">
Labels:
<span class="Ins">Insertion</span>
<span class="Del">Deletion</span>
<span class="Sub">Substitution</span>
<span class="Splice">Splice</span>
<span class="Splice">Inversion</span>
<span class="Unknown">Unknown</span>
<span class="Ins_Allele">Insertion of the predicted allele</span>
<span class="Del_Allele">Deletion of the predicted allele</span>
<span class="Inv_Allele">Inversion of the predicted allele</span>
</p>
<hr>
"""

HTML_FOOTER = """
</body>
</html>
"""


# Build html
def embed_sv_allele_to_html(misdv_sv_allele: list[str]) -> list[str]:
    html_sv_allele = []
    idx = 0
    while idx < len(misdv_sv_allele):
        midsv_tag = misdv_sv_allele[idx]
        if midsv_tag.startswith("+"):
            html_sv_allele.append("<span class='Ins_Allele'>")
            insertions, last_tag = midsv_tag.split("|")[:-1], midsv_tag.split("|")[-1]
            ins_seq = ["=" + ins[1] for ins in insertions]
            html_sv_allele.extend(ins_seq)
            html_sv_allele.append("</span>")
            if last_tag.startswith("-"):
                html_sv_allele.append("<span class='Del_Allele'>")
                html_sv_allele.append(last_tag)
                html_sv_allele.append("</span>")
            else:  # match or substitution
                html_sv_allele.append(last_tag)
        elif midsv_tag.startswith("-"):
            html_sv_allele.append("<span class='Del_Allele'>")
            html_sv_allele.append(midsv_tag)
            html_sv_allele.append("</span>")
        elif midsv_tag.islower():
            html_sv_allele.append("<span class='Inv_Allele'>")
            html_sv_allele.append(midsv_tag)
            html_sv_allele.append("</span>")
        else:
            html_sv_allele.append(midsv_tag)

        idx += 1
    return html_sv_allele


def embed_mutations_to_html(html_sv_allele: list[str], midsv_consensus: list[str]) -> str:
    idx = 0
    for i, tag_sv in enumerate(html_sv_allele):
        if tag_sv.startswith("<"):
            continue
        tag_sample = midsv_consensus[idx]
        if tag_sample.startswith("*"):
            substitutions = [tag_sample[-1]]
            while idx < len(midsv_consensus) - 1 and midsv_consensus[idx + 1].startswith("*"):
                substitutions.append(midsv_consensus[idx + 1][-1])
                idx += 1
            html_sv_allele[i] = "<span class='Sub'>" + "".join(substitutions) + "</span>"
        elif tag_sample.startswith("+"):
            html_sv_allele.append("<span class='Ins'>")
            insertions, last_tag = tag_sample.split("|")[:-1], tag_sample.split("|")[-1]
            ins_seq = "".join(ins[-1] for ins in insertions)
            html_sv_allele.append(ins_seq)
            html_sv_allele.append("</span>")
            if last_tag.startswith("-"):
                html_sv_allele.append("<span class='Del'>")
                html_sv_allele.append(last_tag[-1])
                html_sv_allele.append("</span>")
            else:  # match or substitution
                html_sv_allele.append(last_tag[-1])
        elif tag_sample.startswith("-"):
            html_sv_allele.append("<span class='Del'>")
            html_sv_allele.append(tag_sample[-1])
            html_sv_allele.append("</span>")
        else:
            html_sv_allele[i] = tag_sv[-1]

        idx += 1
    return "".join(html_sv_allele)


def to_html(midsv_sv_allele: list[str], midsv_consensus: list[str:], description: str = "") -> str:
    """Output HTML string showing a sequence with mutations colored"""
    description_str = f"<h1>{description}</h1>" if description else ""
    html_sv_allele = embed_sv_allele_to_html(midsv_sv_allele)
    html_parts = embed_mutations_to_html(html_sv_allele, midsv_consensus)
    return "\n".join(
        [
            HTML_HEADER,
            description_str,
            HTML_LEGEND,
            html_parts,
            HTML_FOOTER,
        ]
    )
