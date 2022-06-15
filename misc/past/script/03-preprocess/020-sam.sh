#!/bin/sh

#----------------------------------------------------------
timestamp "Generate SAM files" >>log_DAJIN.txt
#----------------------------------------------------------

outdir=.DAJIN_temp/sam/
mkdir -p "$outdir" /tmp/"$outdir"/

find .DAJIN_temp/fasta/*.fa |
  while read -r fasta; do
    allele_name="$(basename ${fasta%.fa})"

    if find /tmp/"$outdir"/"$control_name"_"$allele_name".bam 1>/dev/null 2>&1; then
      samtools view /tmp/"$outdir"/"$control_name"_"$allele_name".bam \
        >"$outdir"/"$control_name"_"$allele_name".sam
    else
      minimap2 -t "$threads" --cs=long \
        -ax map-ont "$fasta" .DAJIN_temp/fastq/"$control_name".fq 2>/dev/null \
        >"$outdir"/"$control_name"_"$allele_name".sam
      samtools sort -@ "$threads" "$outdir"/"$control_name"_"$allele_name".sam \
        >/tmp/"$outdir"/"$control_name"_"$allele_name".bam 2>/dev/null
    fi

    minimap2 -t "$threads" --cs=long \
      -ax map-ont "$fasta" .DAJIN_temp/fastq/"$sample_name".fq 2>/dev/null \
      >"$outdir"/"$sample_name"_"$allele_name".sam
  done
