import re
import os
import subprocess


def split_fasta(fasta_file: str, fasta_dir: str) -> None:

    with open(fasta_file, "r") as f:
        regex = re.compile("(>.*?)\n([A-Za-z\n]*)", re.DOTALL)
        fasta_wrap = re.findall(regex, f.read())
        fasta_headers = [f[0].replace(">", "") for f in fasta_wrap]
        fasta_headers = fasta_headers.replace("\t", " ")
        fasta_headers = fasta_headers.replace(",", "_")
        fasta_contents = ["\n".join(f) for f in fasta_wrap]

    for head, content in zip(fasta_headers, fasta_contents):
        output = os.path.join(fasta_dir, head) + ".fasta"
        with open(output, 'w') as f:
            f.write(content)


def minimap2(fastq_file, fasta_dir, sam_dir, threads=1):
    fastq_name = os.path.basename(fastq_file).split(".", 1)[0]
    fasta_files = [os.path.join(fasta_dir, f) for f in os.listdir(fasta_dir)]
    fasta_names = [f.replace(".fasta", "") for f in os.listdir(fasta_dir)]

    output_names = [fastq_name + "_" + f for f in fasta_names]
    for name, fasta in zip(output_names, fasta_files):
        output_sam = os.path.join(sam_dir, name) + ".sam"
        with open(output_sam, "w") as f:
            minimap2 = ["minimap2", "-ax", "map-ont", "--cs=long",
                        "-t", str(threads), fasta, fastq_file]
            subprocess.run(minimap2,
                           stdout=f,
                           stderr=subprocess.DEVNULL,
                           check=True)


# fastq_file = "tests/data/query.fq"
# fasta_dir = os.path.join(".tmpDAJIN", "fasta")
# output_dir = os.path.join(".tmpDAJIN", "sam")
# threads = 10

# minimap2(fastq_file, fasta_dir, output_dir, threads)
