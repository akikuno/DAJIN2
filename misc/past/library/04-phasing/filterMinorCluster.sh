#!/bin/sh

filterMinorCluster() {
  timestamp "Filter minor alleles" log_DAJIN.txt

  cat - |
    tee .DAJIN_temp/clustering/tmp_"$sample_name" |
    cut -d, -f1-2 |
    sort |
    uniq -c |
    awk '{
      sum+=$1
      num[NR]=$1
      allele[NR]=$2
      } END {
      for(i=1; i<=NR; i++)
        print allele[i],num[i],num[i]/sum*100
    }' |
    awk '$3>1 {print $1}' |
    tr "," "_" |
    sort >.DAJIN_temp/clustering/tmp_allele

  cat .DAJIN_temp/clustering/tmp_"$sample_name" |
    sed "s/,/_/" |
    sort |
    join -t, - .DAJIN_temp/clustering/tmp_allele |
    sed "s/_/,/"

  rm .DAJIN_temp/clustering/tmp_*
}
