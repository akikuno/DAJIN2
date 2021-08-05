#!/bin/sh

#----------------------------------------------------------
timestamp "Generate SAM files" >>log_DAJIN.txt
#----------------------------------------------------------

mkdir -p .DAJIN_temp/sam /tmp/DAJIN/sam

find .DAJIN_temp/fasta/*.fa |
  while read -r fasta; do
    allele_name="$(basename ${fasta%.fa})"

    if find /tmp/DAJIN/sam/"$control_name"_"$allele_name".bam 1>/dev/null 2>&1; then
      samtools view /tmp/DAJIN/sam/"$control_name"_"$allele_name".bam \
        >.DAJIN_temp/sam/"$control_name"_"$allele_name".sam
    else
      minimap2 -t "$threads" -ax map-ont "$fasta" .DAJIN_temp/fastq/"$control_name".fq --cs=long 2>/dev/null \
        >.DAJIN_temp/sam/"$control_name"_"$allele_name".sam
      samtools sort -@ "$threads" .DAJIN_temp/sam/"$control_name"_"$allele_name".sam \
        >/tmp/DAJIN/sam/"$control_name"_"$allele_name".bam 2>/dev/null
    fi

    minimap2 -t "$threads" -ax map-ont "$fasta" .DAJIN_temp/fastq/"$sample_name".fq --cs=long 2>/dev/null \
      >.DAJIN_temp/sam/"$sample_name"_"$allele_name".sam
  done
