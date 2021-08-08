#!/bin/sh

. .DAJIN_temp/library/calcFreqMIDS.sh
mkdir -p .DAJIN_temp/consensus /tmp/DAJIN/consensus

# control =================================================

cat .DAJIN_temp/midsmask/barcode32_control.csv |
  cut -d, -f 2- |
  sed "s/^/control,/" |
  sed "s/,[0-9][0-9]*[MDS]/,I/g" |
  calcFreqMIDS |
  sed "s/,/@/" |
  sort -u -t, |
  cat >.DAJIN_temp/consensus/tmp_control_midsfreq.csv

# sample  =================================================

cat .DAJIN_temp/clustering/"$sample_name".csv |
  cut -d, -f1-2 |
  sort -u
.DAJIN_temp/midsmask/"$sample_name"*_control.csv |
