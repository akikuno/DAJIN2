#!/bin/sh

mkdir -p .DAJIN_temp/fasta/

formatFasta "${alleles}" |
  while read -r line; do
    output="$(echo ${line#>} | cut -d " " -f 1)"
    [ _"$output" = "_wt" ] && output="control"
    echo "$line" |
      tr " " "\n" >.DAJIN_temp/fasta/"$output".fa
  done
