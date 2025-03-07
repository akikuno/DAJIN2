def annotate_sv_allele(cons_midsv_tag: list[str], midsv_sv_allele: list[str]) -> list[str]:
    """Annotate the SV allele to consenses midsv tag."""
    cons_midsv_tag_with_sv_tag = []
    idx = 0
    num_deletion = 0
    num_insertion = 0
    while idx < len(midsv_sv_allele) and idx - num_deletion + num_insertion < len(cons_midsv_tag):
        sv_tag = midsv_sv_allele[idx]

        # Deletion
        if sv_tag.startswith("-"):
            cons_midsv_tag_with_sv_tag.append(sv_tag)
            # Enclose consecutive alleles within a single span.
            num_deletion += 1
            idx += 1
            while idx < len(midsv_sv_allele) and midsv_sv_allele[idx].startswith("-"):
                cons_midsv_tag_with_sv_tag.append(midsv_sv_allele[idx])
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
