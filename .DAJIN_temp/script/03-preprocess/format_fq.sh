#!/bin/sh

mkdir -p .DAJIN_temp/fastq/
control_name="$(echo ${control##*/} | sed "s/\..*$//" | tr " " "_")"
sample_name="$(echo ${sample##*/} | sed "s/\..*$//" | tr " " "_")"

fmt_fq "$control" >.DAJIN_temp/fastq/"$control_name".fq
fmt_fq "$sample" >.DAJIN_temp/fastq/"$sample_name".fq
