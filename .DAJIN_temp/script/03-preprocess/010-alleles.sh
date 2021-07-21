#!/bin/sh

. .DAJIN_temp/library/fmt_fa.sh

fmt_fa "${alleles}" |
  while read -r line; do
    output="$(echo ${line#>} | cut -d " " -f 1)"
    echo "$line" |
      tr " " "\n" >.DAJIN_temp/input/"$output".fa
  done
