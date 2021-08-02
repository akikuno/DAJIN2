#!/bin/sh

scoreToScalar() {
  cat "$1" |
    awk -F, 'BEGIN {OFS=","} {
    sum=0
    for(i=2; i<=NF; i++) sum+=$i
    print $1, log(sum/(NF-1))
  }'
}
