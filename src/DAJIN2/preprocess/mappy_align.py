import cstag
import mappy

complement = {"A": "T", "C": "G", "G": "C", "T": "A"}


def to_sam(PATH_REFFA: str, PATH_QUEFQ: str, cslong: bool = True) -> list:
    """Align seqences using mappy and Convert PAF to SAM

    Args:
        PATH_REFFA (str): Path of reference fasta
        PATH_QUEFQ (str): Path of query fasta/fastq
        cslong (bool, optional): long formatted CS tag if True. Defaults to True.

    Returns:
        list: List of SAM
    """
    # SQ header
    SAM = [f"@SQ\tSN:{n}\tLN:{len(s)}" for n, s, _ in mappy.fastx_read(PATH_REFFA)]
    # Mappy
    ref = mappy.Aligner(PATH_REFFA)
    if not ref:
        raise Exception("ERROR: failed to load fasta file")
    for qname, qseq, qual in mappy.fastx_read(PATH_QUEFQ):
        for hit in ref.map(qseq, cs=True):
            # flag
            if hit.is_primary:
                flag = 0 if hit.strand == 1 else 16
            else:
                flag = 2048 if hit.strand == 1 else 2064
            # Append softclips to CIGAR
            cigar = hit.cigar_str
            if hit.q_st > 0:
                softclip = str(hit.q_st) + "S"
                cigar = softclip + cigar if hit.strand == 1 else cigar + softclip
            if len(qseq) - hit.q_en > 0:
                softclip = str(len(qseq) - hit.q_en) + "S"
                cigar = cigar + softclip if hit.strand == 1 else softclip + cigar
            # Revcomp
            if not hit.strand == 1:
                qseq = "".join(complement[nt] for nt in qseq[::-1])
            qseq = qseq.upper()
            # cslong
            cs = "cs:Z:" + hit.cs
            if cslong:
                cs = cstag.lengthen(hit.cs, cigar, qseq)
            # summarize
            alignment = [
                qname,
                flag,
                hit.ctg,
                hit.r_st + 1,
                hit.mapq,
                cigar,
                "*",
                0,
                0,
                qseq,
                qual,
                cs,
            ]
            alignment = [str(a) for a in alignment]
            SAM.append("\t".join(alignment))
    return SAM