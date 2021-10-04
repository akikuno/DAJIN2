import re
import os
import subprocess
fasta_path = "tests/data/ref.fa"
fastq_path = "tests/data/query.fq"
fasta_dir = ".tmpDAJIN/fasta"
sam_dir = ".tmpDAJIN/sam"
threads = 10
os.makedirs(fasta_dir, exist_ok=True)
os.makedirs(sam_dir, exist_ok=True)


with open(fasta_path, "r") as f:
    regex = re.compile("(>.*?)\n([A-Za-z\n]*)", re.DOTALL)
    fasta_wrap = re.findall(regex, f.read())
    fasta_ids = [f[0].replace(">", "") for f in fasta_wrap]
    fasta_contents = ["\n".join(f) for f in fasta_wrap]

for id, content in zip(fasta_ids, fasta_contents):
    with open(f'{fasta_dir}/{id}.fasta', 'w') as f:
        f.write(content)


def ls_fullpath(d):
    return [os.path.join(d, f) for f in os.listdir(d)]


for id, fasta in zip(fasta_ids, ls_fullpath(fasta_dir)):
    with open(f"{sam_dir}/{id}.sam", "w") as f:
        minimap2 = ["minimap2", "-ax", "map-ont", "--cs=long",
                    "-t", str(threads), fasta, fastq_path]
        subprocess.run(minimap2, stdout=f, check=True,)
