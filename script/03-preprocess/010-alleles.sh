#!/bin/sh

cat "${alleles}" |
  awk 'BEGIN{RS=">"} {for(i=1;i<=NF;i++) print $i}'
sed "s/ /_/g" |
  awk '$1~/^>/ {$1=$1"@SEP@"}1' |
  tr -d "\r\n" |
  sed "s/@SEP@/ /g" |
  awk '{print $1, toupper($2)}' |
  grep -v "^>$" |
  while read -r line; do
    output="$(echo ${line#>} | cut -d " " -f 1)"
    echo "$line" |
      tr " " "\n" >.DAJIN_temp/input/"$output".fa
  done
