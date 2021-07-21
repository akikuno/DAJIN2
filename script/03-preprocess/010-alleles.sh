#!/bin/sh

. .DAJIN_temp/library/fmt_fa.sh

mkdir -p .DAJIN_temp/fasta/"$output".fa

fmt_fa "${alleles}" |
  while read -r line; do
    output="$(echo ${line#>} | cut -d " " -f 1)"
    echo "$line" |
      tr " " "\n" >.DAJIN_temp/input/"$output".fa
  done
