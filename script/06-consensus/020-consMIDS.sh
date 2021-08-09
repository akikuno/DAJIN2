#!/bin/sh

outdir=.DAJIN_temp/consensus/consmids/
mkdir -p "$outdir" /tmp/"$outdir"

#? MEMO =========================
#?  =========================

find .DAJIN_temp/consensus/freqmids/"$sample_name"* |
  while read -r line; do
    outfile="$(basename "$line")"
    cat "$line" |
      join -t, - .DAJIN_temp/consensus/freqmids/"$control_name".csv |
      awk -F, 'function RELU(v) {return v < 0 ? 0 : v}
        BEGIN {OFS=","} {
        subI=RELU($3-$7)
        subD=RELU($4-$8)
        subS=RELU($5-$9)
        subM=RELU(1-subI-subD-subS)
        # print subM,subI,subD,subS
        if ($3>0.9) MIDS="I"
        else if ($4>0.9) MIDS="D"
        else if ($5>0.9) MIDS="S"
        #* sequence error
        else if (subI>subM && subI>subD && subI>subS) MIDS="I"
        else if (subD>subM && subD>subI && subD>subS) MIDS="D"
        else if (subS>subM && subS>subI && subS>subD) MIDS="S"
        else MIDS="M"
        print $0,MIDS
      }' |
      tr "@" "," |
      awk -F, '{print $2","$NF}' |
      sort -t, -n |
      grep ^ >"$outdir"/"$outfile"
  done
