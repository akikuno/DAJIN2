#!/bin/sh

fmt_fa() {
  cat "$1" |
    awk 'BEGIN{RS=">"} {
      $1=">" $1 " "
      for(i=1;i<=NF;i++) printf $i
      print ""
      }' |
    awk 'NF==2'
}
