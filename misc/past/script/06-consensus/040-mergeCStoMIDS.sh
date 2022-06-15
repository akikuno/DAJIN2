#!/bin/sh

outdir=.DAJIN_temp/consensus/merge/
mkdir -p "$outdir"

find .DAJIN_temp/consensus/fmtcs/*.csv |
  while read -r line; do
    consMIDS=.DAJIN_temp/consensus/consmids/"$(basename "$line")"
    if grep -q "[IDS]" "$consMIDS"; then
      Rscript --vanilla "$(find .DAJIN_temp -name margeCStoMIDS.R)" \
        "$line" "$consMIDS" |
        cat >"$outdir"/"$(basename "$line")"
    fi
  done

#? DEBUG ============================================
# cat "$line" >tmp_cs
# cat "$consMIDS" >tmp_mids
