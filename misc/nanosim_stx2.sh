#!/bin/bash

mkdir -p examples/nanosim/del-stx2
mkdir -p examples/nanosim/inv-stx2

design="examples/del-stx2/design_stx2.fa"
cat "$design" |
    paste - - |
    grep control |
    tr "\t" "\n" |
    cat >tmp_ref.fa

cat "$design" |
    paste - - |
    grep target |
    tr "\t" "\n" |
    cat >tmp_del.fa

cat "$design" |
    paste - - |
    grep inversion |
    tr "\t" "\n" |
    cat >tmp_inv.fa

control=examples/del-stx2/barcode30.fq.gz
threads=20
conda activate nanosim

read_analysis.py genome -t "$threads" -i "$control" -rg tmp_ref.fa

simulator.py genome -t "$threads" -n 3000 -rg tmp_ref.fa -b guppy --fastq -o control
simulator.py genome -t "$threads" -n 3000 -rg tmp_del.fa -b guppy --fastq -o del
simulator.py genome -t "$threads" -n 3000 -rg tmp_inv.fa -b guppy --fastq -o inv

# リード数の確認
for fq in *_aligned_reads.fastq; do
    echo "$fq"
    grep -c "^@" "$fq"
done

# IGVで可視化

for fq in *_aligned_reads.fastq; do
    echo "$fq"
    minimap2 -ax map-ont -t "$threads" tmp_ref.fa "$fq" |
        samtools sort -@ "$threads" >"${fq%.fastq}".bam
    samtools index -@ "$threads" "${fq%.fastq}".bam
done

# 確認後、fastqを保存
gzip -c control_aligned_reads.fastq >examples/nanosim/del-stx2/control.fq.gz
gzip -c del_aligned_reads.fastq >examples/nanosim/del-stx2/deletion.fq.gz

gzip -c control_aligned_reads.fastq >examples/nanosim/inv-stx2/control.fq.gz
gzip -c inv_aligned_reads.fastq >examples/nanosim/inv-stx2/inversion.fq.gz

rm tmp* training* simulated* *aligned_reads* *error*
