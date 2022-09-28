#!/bin/bash

mkdir -p examples/nanosim/pm-tyr
cp examples/pm-tyr/design_tyr.fa examples/nanosim/pm-tyr/design_tyr.fa

cat examples/pm-tyr/design_tyr.fa |
    paste - - |
    grep control |
    tr "\t" "\n" |
    cat >tmp_ref.fa

cat examples/pm-tyr/design_tyr.fa |
    paste - - |
    grep albino |
    tr "\t" "\n" |
    cat >tmp_query.fa

control=examples/pm-tyr/barcode32.fq.gz
threads=20

conda activate nanosim

read_analysis.py genome -t "$threads" -i "$control" -rg tmp_ref.fa

simulator.py genome -t "$threads" -n 1000 -rg tmp_ref.fa -b guppy --fastq -o control
simulator.py genome -t "$threads" -n 1000 -rg tmp_query.fa -b guppy --fastq -o albino

gzip -c control_aligned_reads.fastq >examples/nanosim/pm-tyr/control.fq.gz
gzip -c albino_aligned_reads.fastq >examples/nanosim/pm-tyr/albino.fq.gz

# minimap2 -ax map-ont tmp_ref.fa control_aligned_reads.fastq |
#     samtools sort >tmp_ref.bam
# samtools index tmp_ref.bam

# minimap2 -ax map-ont tmp_ref.fa albino_aligned_reads.fastq |
#     samtools sort >tmp_query.bam
# samtools index tmp_query.bam

rm tmp* training* simulated* albino* control*
