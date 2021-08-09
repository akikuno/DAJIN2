#!/bin/sh

. "$(find .DAJIN_temp -name calcFreqMIDS.sh)"

outdir=.DAJIN_temp/consensus/freqcs/
mkdir -p "$outdir" /tmp/"$outdir"

#? MEMO =========================
#? calcFreqMIDSの１列目のID指定はいらないみたい
#?  =========================

find .DAJIN_temp/consensus/fmtcs/*.csv |
  while read -r line; do
    output="$outdir"/"$(basename "$line")"
    consmids=.DAJIN_temp/consensus/consmids/"$(basename "$line")"

    cat "$line" |
      cut -d, -f 2- |
      cat >tmp_cs

    cat $consmids >tmp_mids
    # calcFreqMIDS |
    #   sed "s/,/@/" |
    #   sort -u -t, |
    #   cat >"$outdir"/"$sample_name"_"$suffix".csv
  done

rm "$outdir"/*
