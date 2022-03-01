import os
import sys
from src.DAJIN2.utils import io

########################################################################
# アレルを辞書型およびシングルfastaフォーマットに変換
########################################################################


def save_allele_as_single_fasta_files(allele: str) -> None:
    fasta = io.fread(allele)
    header = []
    sequence = []
    idx = -1
    for f in fasta:
        if ">" in f:
            header.append(f)
            sequence.append([])
            idx += 1
        else:
            sequence[idx] += f
    sequence = ["".join(s).upper() for s in sequence]
    if len(header) > len(set(header)):
        print(
            f"Error: '{allele}' includes duplicated sequences.\n"
            f"'{allele}' must contain only unique DNA sequences.",
            file=sys.stderr,
        )
        sys.exit(1)
    if header.count(">control") == 0:
        print(
            f"Error: '{allele}' does not contain 'control' sequence.\n"
            f"'{allele}' must include a 'control' sequence.",
            file=sys.stderr,
        )
        sys.exit(1)
    for h, s in zip(header, sequence):
        filename = h.replace(">", "")
        contents = "\n".join([h, s]) + "\n"
        with open(f".tmpDAJIN/fasta/{filename}.fasta", "w") as f:
            f.write(contents)


def dictionize_allele(allele: str) -> dict:
    fasta = io.fread(allele)
    header = []
    sequence = []
    idx = -1
    for f in fasta:
        if ">" in f:
            header.append(f)
            sequence.append([])
            idx += 1
        else:
            sequence[idx] += f
    sequence = ["".join(s).upper() for s in sequence]
    return {h.replace(">", ""): s for h, s in zip(header, sequence)}


dict_allele = dictionize_allele(allele)

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


for fastqpath in [sample, control]:
    fastq_anno = annotate_TooLong_to_fastq(fastqpath, dict_allele)
    basename = os.path.basename(fastqpath)
    io.fwrite(fastq_anno, f".tmpDAJIN/fastq/{basename}")

