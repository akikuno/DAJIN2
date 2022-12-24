from copy import deepcopy
import mappy
import cstag
from pathlib import Path

path_reference_fasta = "tests/data/mappy/stx2_control.fa"
path_query_fastq = "tests/data/mappy/stx2_del.fq"
cslong = True


def revcomp(sequence: str) -> str:
    complement = {"A": "T", "C": "G", "G": "C", "T": "A"}
    return "".join(complement[nt] for nt in sequence[::-1])


SAM = [f"@SQ\tSN:{n}\tLN:{len(s)}" for n, s, _ in mappy.fastx_read(path_reference_fasta)]
ref = mappy.Aligner(path_reference_fasta)

if not ref:
    raise AttributeError(f"Failed to load f{path_reference_fasta}")

for MAPPY_NAME, MAPPY_SEQ, MAPPY_QUAL in mappy.fastx_read(path_query_fastq):
    for hit in ref.map(MAPPY_SEQ, cs=True):
        query_seq = deepcopy(MAPPY_SEQ)
        query_qual = deepcopy(MAPPY_QUAL)
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
        if len(MAPPY_SEQ) - hit.q_en > 0:
            softclip = str(len(MAPPY_SEQ) - hit.q_en) + "S"
            cigar = cigar + softclip if hit.strand == 1 else softclip + cigar
        # Revcomp
        if hit.strand == -1:
            query_seq = revcomp(query_seq)
            query_qual = query_qual[::-1]
        query_seq = query_seq.upper()
        # cslong
        cs = "cs:Z:" + hit.cs
        if cslong:
            cs = cstag.lengthen(hit.cs, cigar, query_seq)
        # summarize
        alignment = [
            MAPPY_NAME,
            flag,
            hit.ctg,
            hit.r_st + 1,
            hit.mapq,
            cigar,
            "*",
            0,
            0,
            query_seq,
            query_qual,
            cs,
        ]
        alignment = [str(a) for a in alignment]
        SAM.append("\t".join(alignment))


Path("tmp.sam").write_text("\n".join(SAM))
