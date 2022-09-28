#!/bin/bash

conda activate DAJIN2

cat examples/del-stx2/design_stx2.fa |
    grep "control" -A1 >tmp_ref.fa

find examples/del-stx2/barcode*.fq.gz |
    while read -r fq; do
        output="$(basename ${fq%.fq.gz})"
        minimap2 -t "${threads:-1}" -ax map-ont tmp_ref.fa "$fq" |
            samtools sort -@ "${threads:-1}" >tmp_"$output".bam
        samtools index -@ "${threads:-1}" tmp_"$output".bam
    done
