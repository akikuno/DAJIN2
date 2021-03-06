#!/bin/bash

# cat barcode41.fastq | head -n 5000 | gzip -c >barcode41.fq.gz
# cat barcode47.fastq | head -n 5000 | gzip -c >barcode47.fq.gz
# cat barcode48.fastq | head -n 5000 | gzip -c >barcode48.fq.gz

directory="examples/flox-cables2/AyabeTask2"
chromosome="chr2"
threads=30

wget -O - https://hgdownload.soe.ucsc.edu/goldenPath/mm39/chromosomes/"$chromosome".fa.gz |
    gzip -dc >tmp_chr.fa

mkdir -p "$directory"/mm39

find "$directory"/barcode*.fq.gz |
    while read -r fq; do
        output="${fq%.fq.gz}".bam
        minimap2 -t "$threads" -ax map-ont tmp_chr.fa "$fq" --cs=long |
            samtools sort -@ "$threads" >"$output"
        samtools index -@ "$threads" "$output"
        mv "$output"* "$directory"/mm39
    done

rm tmp_*
