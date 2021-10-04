import re
import os
import subprocess


def split_fasta(fasta_file: str, output_dir: str) -> None:

    with open(fasta_file, "r") as f:
        regex = re.compile("(>.*?)\n([A-Za-z\n]*)", re.DOTALL)
        fasta_wrap = re.findall(regex, f.read())
        fasta_headers = [f[0].replace(">", "") for f in fasta_wrap]
        fasta_contents = ["\n".join(f) for f in fasta_wrap]

    for head, content in zip(fasta_headers, fasta_contents):
        output = os.path.join(output_dir, head) + ".fasta"
        with open(output, 'w') as f:
            f.write(content)


def minimap2(fasta_dir, fastq_file, output_dir, threads=1):
    fasta_files = [os.path.join(fasta_dir, f) for f in os.listdir(fasta_dir)]
    heads = [f.replace(".fasta", "") for f in os.listdir(fasta_dir)]
    for head, fasta in zip(heads, fasta_files):
        output_sam = os.path.join(output_dir, head) + ".sam"
        with open(output_sam, "w") as f:
            minimap2 = ["minimap2", "-ax", "map-ont", "--cs=long",
                        "-t", str(threads), fasta, fastq_file]
            subprocess.run(minimap2, stdout=f, check=True)


# fastq_file = "tests/data/query.fq"
# fasta_dir = os.path.join(".tmpDAJIN", "fasta")
# output_dir = os.path.join(".tmpDAJIN", "sam")
# threads = 10

# minimap2(fasta_dir, fastq_file, output_dir, threads)
