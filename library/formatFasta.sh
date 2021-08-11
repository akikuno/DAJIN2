#!/bin/sh

formatFasta() {
  cat "$1" |
    tr -d "\r" |
    awk 'BEGIN {RS = ">"} {
      $1=">" $1 " "
      for(i=1;i<=NF;i++) printf $i
      print ""
      }' |
    awk 'NF==2 {$2=toupper($2)}1'
}
