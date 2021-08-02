#!/bin/sh

. .DAJIN_temp/library/formatFasta.sh

mkdir -p .DAJIN_temp/fasta/

formatFasta "${alleles}" |
  while read -r line; do
    output="$(echo ${line#>} | cut -d " " -f 1)"
    echo "$line" | tr " " "\n" >.DAJIN_temp/fasta/"$output".fa
  done
