#!/bin/bash

threads=30

directory="examples/del-stx2"
chromosome="chr5"

wget -O - https://hgdownload.soe.ucsc.edu/goldenPath/mm39/chromosomes/"$chromosome".fa.gz |
    gzip -dc >tmp_chr.fa

mkdir -p "$directory"/mm39

find "$directory"/barcode*.fq.gz |
    while read -r fq; do
        output="${fq%.fq.gz}".bam
        minimap2 -t "$threads" -ax map-ont tmp_chr.fa "$fq" --cs=long |
            samtools sort -@ "$threads" >"$output"
        samtools index "$output"
        mv "$output"* "$directory"/mm39
    done

rm tmp_*
