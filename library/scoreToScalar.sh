#!/bin/sh

scoreToScalar() {
  cat "$1" |
    awk -F, '{
    sum=0
    for(i=2; i<=NF; i++) {
      sum+=$i
    }
    print log(sum/(NF-1))
  }'
}
