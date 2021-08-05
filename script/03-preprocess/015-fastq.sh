#!/bin/sh

mkdir -p .DAJIN_temp/fastq/ /tmp/DAJIN/fastq/

# if find /tmp/fastq/"$control_name"* 1>/dev/null 2>&1; then
#   gzip -dc /tmp/fastq/"$control_name".fq.gz >.DAJIN_temp/fastq/"$control_name".fq
#   open "$sample" | maskFastq >.DAJIN_temp/fastq/"$sample_name".fq
# else
#   {
#     open "$control" | maskFastq >.DAJIN_temp/fastq/"$control_name".fq &
#     open "$sample" | maskFastq >.DAJIN_temp/fastq/"$sample_name".fq &
#     wait
#   } 1>/dev/null 2>&1
#   gzip -c .DAJIN_temp/fastq/"$control_name".fq >/tmp/fastq/"$control_name".fq.gz
# fi

open "$control" >.DAJIN_temp/fastq/"$control_name".fq
open "$sample" >.DAJIN_temp/fastq/"$sample_name".fq
