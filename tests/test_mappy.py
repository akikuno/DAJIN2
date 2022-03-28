import mappy as mp

ref_fasta = "tests/data/mappy/ref.fa"
que = "tests/data/mappy/query.fq"

ref_name, ref_seq, _ = list(mp.fastx_read(ref_fasta))[0]

ref = mp.Aligner(ref_fasta)

tmp = []
for name, seq, qual in mp.fastx_read(que):
    for hit in ref.map(seq, cs=True):
        tmp.append([hit.r_st, hit.q_st, hit.cigar_str, hit.cs])

tmp
for t in tmp:
    if "S" in t[1]:
        break
