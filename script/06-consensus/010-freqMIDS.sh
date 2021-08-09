#!/bin/sh

. .DAJIN_temp/library/calcFreqMIDS.sh

outdir=.DAJIN_temp/consensus/freqmids/
mkdir -p "$outdir" /tmp/"$outdir"

#? MEMO =========================
#? calcFreqMIDSの１列目のID指定はいらないみたい
#?  =========================

# control =================================================

if [ -r /tmp/"$outdir"/"$control_name".csv.gz ]; then
  cp /tmp/"$outdir"/"$control_name".csv.gz "$outdir"
  gzip -f -d "$outdir"/"$control_name".csv.gz
else
  cat .DAJIN_temp/midsmask/"$control_name"_control.csv |
    cut -d, -f 2- |
    sed "s/^/control,/" |
    sed "s/,[0-9][0-9]*[MDS]/,I/g" |
    calcFreqMIDS |
    sed "s/,/@/" |
    sort -u -t, |
    cat >"$outdir"/"$control_name".csv
  cp "$outdir"/"$control_name".csv /tmp/"$outdir"/
  gzip /tmp/"$outdir"/"$control_name".csv
fi

# sample  =================================================

cat .DAJIN_temp/clustering/"$sample_name".csv |
  cut -d, -f1-2 |
  sort -u |
  while read -r line; do
    suffix="$(echo "$line" | tr , _)"
    grep "^$line" .DAJIN_temp/clustering/"$sample_name".csv |
      cut -d, -f3 |
      sort -u >.DAJIN_temp/clustering/tmp_"$suffix"

    cat .DAJIN_temp/midsmask/"$sample_name"*_control.csv |
      join -t, - .DAJIN_temp/clustering/tmp_"$suffix" |
      cut -d, -f 2- |
      sed "s/^/control,/" |
      sed "s/,[0-9][0-9]*[MDS]/,I/g" |
      calcFreqMIDS |
      sed "s/,/@/" |
      sort -u -t, |
      cat >"$outdir"/"$sample_name"_"$suffix".csv
  done
