import mappy


def mappy2sam(REFFA: str, QUEFQ: str) -> list:
    # SQ header
    SAM = [f"@SQ\tSN:{n}\tLN:{len(s)}" for n, s, _ in mappy.fastx_read(REFFA)]
    # Mappy
    ref = mappy.Aligner(REFFA)
    if not ref:
        raise Exception("ERROR: failed to load fasta file")
    for qname, qseq, qual in mappy.fastx_read(QUEFQ):
        for hit in ref.map(qseq, cs=True):
            # flag
            if hit.is_primary:
                flag = 0 if hit.strand == 1 else 16
            else:
                flag = 2048 if hit.strand == 1 else 2064
            # Append softclips to CIGAR
            cigar = hit.cigar_str
            if hit.q_st > 0:
                cigar = str(hit.q_st) + "S" + cigar
            if len(qseq) - hit.q_en > 0:
                cigar = cigar + str(len(qseq) - hit.q_en) + "S"
            # summarize
            alignment = [qname, flag, hit.ctg, hit.r_st + 1, hit.mapq, cigar, "*", 0, 0, qseq, qual, "cs:Z:" + hit.cs]
            alignment = [str(a) for a in alignment]
            SAM.append("\t".join(alignment))
    return SAM
