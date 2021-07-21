#!/bin/sh

cat "${alleles}" |
  sed "s/ /_/g" |
  tr -d "\r" |
  paste - - |
  awk '{print $1, toupper($2)}' |
  grep -v "^>$" |
  while read -r line; do
    output="$(echo ${line#>} | cut -d " " -f 1)"
    echo "$line" |
      tr " " "\n" >.DAJIN_temp/input/"$output".fa
  done
