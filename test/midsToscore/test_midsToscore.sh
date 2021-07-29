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
bbb,M,2M,D,D,M
ccc,M,1M,M,M,S
EOF
  expansion >tmp.csv

# cat <<EOF |
# aaa,M,2M,M,D,M,M,M,M,M
# bbb,M,2M,D,D,M,D,D,D,M
# ccc,M,1M,M,D,S,S,M,M,M
# EOF
#   expansion >tmp.csv

# cat <<EOF |
# aaa,D,D,M,D,D,M
# bbb,D,D,M,D,D,M
# EOF
#   expansion >tmp.csv

cat tmp.csv |
  awk -f library/mutToat.awk |
  awk -f library/atTonum.awk >tmp_row.csv

cat tmp.csv |
  Rscript library/colTable.R >tmp_col.csv

cat tmp.csv |
  awk '{
    FS=","
    OFS=""
    $1=$1","
    for(i=2;i<=NF;i++) {
      if($i>0 && $(i+1)>0)
        $i="@"
      else if($(i-1)=="@" && $i==1)
        $i="@,"
      else
        $i=$i","
    }
    sub(",$", "", $NF)
  }1'

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
