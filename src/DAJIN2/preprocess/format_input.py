# sample, control, allele, output, genome, debug, threads = (
#     "examples/pm-tyr/barcode31.fq.gz",
#     "examples/pm-tyr/barcode32.fq.gz",
#     "examples/pm-tyr/design_tyr.fa",
#     "DAJIN_results",
#     "mm10",
#     True,
#     14
#     )

import gzip
import os
import sys


def fread(filepath: str) -> list:
    """
    read gzip or text file
    """
    with gzip.open(filepath, "rt") as f:
        try:
            out = f.read().splitlines()
        except:
            with open(filepath) as f:
                out = f.read().splitlines()
    return out


def fwrite(infile: list, outfilepath: str) -> None:
    """
    save gzipped file
    """
    with gzip.open(outfilepath, "wt") as f:
        for i in infile:
            f.write("\n".join(i) + "\n")


########################################################################
# アレルを辞書型およびシングルfastaフォーマットに変換
########################################################################


def dictionize_allele(allele: str) -> dict:
    fasta = fread(allele)
    key = []
    seq = []
    s = []
    for f in fasta:
        if ">" in f:
            key.append(f[1:])
            seq.append("".join(s))
        else:
            s.append(f.upper())
    seq.append("".join(s))
    seq = seq[1:]
    if len(key) > len(set(key)):
        print(
            f"Error: '{allele}' includes duplicated sequences.\n'{allele}' must contain only unique DNA sequences.",
            file=sys.stderr,
        )
        sys.exit(1)
    if key.count("control") == 0:
        print(
            f"Error: '{allele}' does not contain 'control' sequences.\n'{allele}' must include 'control' or 'wt' sequences.",
            file=sys.stderr,
        )
        sys.exit(1)
    return {k: s for k, s in zip(key, seq)}


dict_allele = dictionize_allele(allele)
control_length = len(dict_allele["hoge"])

########################################################################
# Fastqのうち過剰に長い配列にTooLongフラグを立てる
########################################################################


for file in [sample, control]:
    fastq = fread(file)
    fastq = [fastq[i : i + 4] for i in range(0, len(fastq), 4)]
    fastq_anno = []
    for f in fastq:
        header = f[0].split()[0]
        if len(f[1]) > control_length * 1.1:
            header = "-".join([header, "TooLong"])
        fastq_anno.append([header, *f[1:]])
    p = os.path.basename(file)
    fwrite(fastq_anno, f".tmpDAJIN/fastq/{p}")

