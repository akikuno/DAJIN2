#!/bin/sh

mkdir -p .DAJIN_temp/midsfilter/

find .DAJIN_temp/mids/*.csv |
  while read -r line; do
    cat "$line" |
      awk -F, '{
      leftD=0
      rightD=0
      for (i=2; i<=51; i++) {
        if ($i=="D") leftD++
      }
      for (i=NF-49; i<=NF; i++) {
        if ($i=="D") rightD++
      }
      if ( leftD!=50 && rightD!=50)
        print $0
      }' >.DAJIN_temp/midsfilter/"$(basename "$line")"
  done
