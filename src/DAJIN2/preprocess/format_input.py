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


def format_allele(allele: str) -> dict:
    """
    - Dictionize allele fasta file
    - Save allele fasta file at ".tmpDAJIN/fasta" directory
    """
    fasta = fread(allele)
    key = []
    seq = []
    s = []
    for f in fasta:
        if ">" in f:
            key.append(f)
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
    if key.count(">control") == 0:
        print(
            f"Error: '{allele}' does not contain 'control' sequences.\n'{allele}' must include 'control' or 'wt' sequences.",
            file=sys.stderr,
        )
        sys.exit(1)
    for k, s in zip(key, seq):
        with open(f".tmpDAJIN/fasta/{k[1:]}.fa", "w") as f:
            f.write("\n".join([k, s]) + "\n")
    return {k[1:]: s for k, s in zip(key, seq)}


dict_allele = format_allele(allele)

########################################################################
# Fastqのうち過剰に長い配列にTooLongフラグを立てる
########################################################################


def format_fastq(fastqpath: str, dict_allele: dict) -> None:
    """
    - Annotate "TooLong" tag after a header when sequence lenth is more than 1.1 control length
    - Save fastq file at ".tmpDAJIN/fastq" directory
    """
    fastq = fread(fastqpath)
    control_length = len(dict_allele["control"])
    # Combine four lines into one list
    fastq = [fastq[i : i + 4] for i in range(0, len(fastq), 4)]
    fastq_anno = []
    for f in fastq:
        header = f[0].split()[0]
        if len(f[1]) > control_length * 1.1:
            header = "-".join([header, "TooLong"])
        fastq_anno.append([header, *f[1:]])
    basename = os.path.basename(fastqpath)
    fwrite(fastq_anno, f".tmpDAJIN/fastq/{basename}")

