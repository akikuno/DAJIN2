#!/bin/sh

#----------------------------------------------------------
timestamp "Generate BAM files" log_DAJIN.txt
#----------------------------------------------------------

mkdir -p .DAJIN_temp/bam/ /tmp/DAJIN/bam/

# Control
if find /tmp/DAJIN/bam/"$control_name".bam >/dev/null 2>&1; then
  cp /tmp/DAJIN/bam/"$control_name".bam* .DAJIN_temp/bam/
else
  cat .DAJIN_temp/sam/"$control_name"*_control.sam |
    samtools sort -@ "$threads" >.DAJIN_temp/bam/"$control_name".bam 2>/dev/null
  samtools index -@ "$threads" .DAJIN_temp/bam/"$control_name".bam
  cp .DAJIN_temp/bam/"$control_name".bam /tmp/DAJIN/bam/
fi

# Sample
cat .DAJIN_temp/sam/"$sample_name"*_control.sam |
  samtools sort -@ "$threads" >.DAJIN_temp/bam/"$sample_name".bam 2>/dev/null
samtools index -@ "$threads" .DAJIN_temp/bam/"$sample_name".bam

cat .DAJIN_temp/clustering/"$sample_name".csv |
  cut -d, -f 1,2 |
  sort -u |
  while read -r line; do
    classif="${line%%,*}"
    number="${line##*,}"
    suffix="$classif"_"$number"

    grep "^$line" .DAJIN_temp/clustering/"$sample_name".csv |
      cut -d, -f3 |
      sort -u >.DAJIN_temp/bam/tmp_id

    grep "^@" .DAJIN_temp/sam/"$sample_name"*_control.sam >.DAJIN_temp/bam/tmp_header

    grep -v "^@" .DAJIN_temp/sam/"$sample_name"*_control.sam |
      sort |
      join - .DAJIN_temp/bam/tmp_id |
      tr " " "\t" |
      cat .DAJIN_temp/bam/tmp_header - |
      samtools sort -@ "$threads" >.DAJIN_temp/bam/"$sample_name"_"$suffix".bam 2>/dev/null
    samtools index -@ "$threads" .DAJIN_temp/bam/"$sample_name"_"$suffix".bam
  done

rm .DAJIN_temp/bam/tmp_*
