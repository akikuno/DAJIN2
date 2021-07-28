#!/bin/sh

. library/midsToscore.sh

echo "aaa,M,M,M,M" |
  midsToscore

echo "aaa,S,S,S,S,S" |
  midsToscore

echo "aaa,M,100S,D,M" |
  midsToscore

echo "aaa,D,100S,D,D" |
  midsToscore

cat <<EOF |
aaa,M,2M,M,D,M
bbb,M,2M,D,D,D
ccc,M,1M,M,D,S
EOF
  expansion >tmp.csv

rowSums tmp

cat tmp |
  awk -F, '{
    id[NR]=$1
    for(i=2; i<=NF; i++) {
      sums[i]+=$i
      if($i>0) {
        array[i, $i]++
      } else {
        array[i, $i]=0
      }
    }
  } END {
    for(nr=1; nr<=NR; nr++) {
      printf id[nr]","
      for(i=2; i<=NF; i++) {
        for(sum=1; sum<=sums[i]; sum++) {
          printf array[i, sum]
        }
        }
      }
    }'
