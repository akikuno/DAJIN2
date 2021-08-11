#!/bin/sh

scoreToScalar() {
  if [ -p /dev/stdin ] && [ "$#" -eq 0 ]; then
    cat -
  elif [ -r "$1" ]; then
    cat "$1"
  else
    echo "$*"
  fi |
    awk -F, 'BEGIN {OFS=","} {
    sum=0
    for(i=2; i<=NF; i++) sum+=$i
    print $1, log(sum/(NF-1))
  }'
}
