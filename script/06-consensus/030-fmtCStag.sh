#!/bin/sh

outdir=.DAJIN_temp/consensus/conscs/
mkdir -p "$outdir" /tmp/"$outdir"

#? MEMO =========================
#? CSタグのコンセンサス配列をつくります
#?  =========================

cat .DAJIN_temp/sam/"$sample_name"_control.sam |
  head
