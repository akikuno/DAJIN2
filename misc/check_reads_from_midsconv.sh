#!/bin/bash

sample=barcode25
readid=f44a2a8c-5a9a-4b94-86da-7f19529cc8df

find .tmpDAJIN/midsconv/"$sample"* |
    while read -r line; do
        echo $line
        grep "$readid" "$line"
    done

cat .tmpDAJIN/sam/barcode25_target.sam |
    grep -e "^@" -e "$readid" |
    samtools sort >tmp.bam
samtools index tmp.bam

cat .tmpDAJIN/sam/barcode25_control.sam |
    grep -e "^@" -e "$readid" |
    samtools sort >tmp2.bam
samtools index tmp2.bam

# extract_full_length_readsが仕事していない？
