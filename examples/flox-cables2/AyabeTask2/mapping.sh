#!/bin/bash

# cat barcode41.fastq | head -n 5000 | gzip -c >barcode41.fq.gz
# cat barcode47.fastq | head -n 5000 | gzip -c >barcode47.fq.gz
# cat barcode48.fastq | head -n 5000 | gzip -c >barcode48.fq.gz

threads=30

wget -O - https://hgdownload.soe.ucsc.edu/goldenPath/mm39/chromosomes/chr2.fa.gz |
    gzip -dc >tmp_chr.fa

mkdir -p mm39

find barcode*.fq.gz |
    while read -r fq; do
        output="${fq%.fq.gz}".bam
        minimap2 -t "$threads" -ax map-ont tmp_chr.fa "$fq" |
            samtools sort -@ "$threads" >mm39/"$output"
        samtools index -@ "$threads" mm39/"$output"
    done

rm tmp_*
