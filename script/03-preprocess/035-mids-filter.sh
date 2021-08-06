#!/bin/sh

mkdir -p .DAJIN_temp/midsfilter/

find .DAJIN_temp/mids/*.csv |
  while read -r line; do
    awk -F, '$51 != "D" && $(NF-50) != "D"' "$line" |
      cat >.DAJIN_temp/midsfilter/"$(basename "$line")"
  done
