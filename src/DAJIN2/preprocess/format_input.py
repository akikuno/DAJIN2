import os
import sys
import mappy
from src.DAJIN2.utils import io

########################################################################
# サンプルが適切なフォーマットなのかチェックする
########################################################################

list(mappy.fastx_read("test.fa"))
list(mappy.fastx_read("test.hogehoge"))

########################################################################
# アレルを辞書型およびシングルfastaフォーマットに変換
########################################################################


def dictionize_allele(allele: str) -> dict:
    header, sequence = [], []
    for name, seq, _ in mappy.fastx_read(allele):
        header.append(name)
        sequence.append(seq.upper())
    # Error handling ------------
    if len(header) > len(set(header)):
        print(
            f"Error: '{allele}' includes duplicated sequences.\n" f"'{allele}' must contain only unique DNA sequences.",
            file=sys.stderr,
        )
        sys.exit(1)
    if header.count("control") == 0:
        print(
            f"Error: '{allele}' does not contain 'control' sequence.\n"
            f"'{allele}' must include a 'control' sequence.",
            file=sys.stderr,
        )
        sys.exit(1)
    # Error handling ------------
    return {h: s for h, s in zip(header, sequence)}


########################################################################
# Fastqのうち過剰に長い配列にTooLongフラグを立てる
########################################################################


def annotate_TooLong_to_fastq(fastqpath: str, dict_allele: dict) -> None:
    """
    Annotate "TooLong" tag after a header when sequence lenth is more than 1.2 times of max sequence length
    """
    fastq = io.fread(fastqpath)
    max_length = max(len(seq) for seq in dict_allele.values())
    # Combine four lines into one list
    fastq = [fastq[i : i + 4] for i in range(0, len(fastq), 4)]
    fastq_anno = []
    for f in fastq:
        header = f[0].split()[0]
        if len(f[1]) > max_length * 1.2:
            header = "-".join([header, "TooLong"])
        fastq_anno.append([header, *f[1:]])
    return fastq_anno
