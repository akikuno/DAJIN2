#!/bin/sh

mkdir -p .DAJIN_temp/midsmaskBySubtract
midsmaskBySubtract="$(find .DAJIN_temp -name "maskBySubtract.R")"

find .DAJIN_temp/midsmaskByPhred/"$sample_name"* >.DAJIN_temp/midsmaskBySubtract/tmp_sample
find .DAJIN_temp/midsmaskByPhred/"$control_name"* >.DAJIN_temp/midsmaskBySubtract/tmp_control

paste .DAJIN_temp/midsmaskBySubtract/tmp_sample .DAJIN_temp/midsmaskBySubtract/tmp_control |
  while read -r line; do
    set $line
    Rscript --vanilla "$midsmaskBySubtract" "$1" "$2" \
      >.DAJIN_temp/midsmaskBySubtract/$(basename "$1")
  done

rm .DAJIN_temp/midsmaskBySubtract/tmp_*

# cp $1 tmp_query.csv
# cp $2 tmp_control.csv
# set --
# cp .DAJIN_temp/midsmask/barcode32_control.csv tmp_c
# Rscript --vanilla --slave "$maskDel" tmp_s tmp_c >tmp_del
# head -n 1 .DAJIN_temp/midsmask/barcode31_control.csv
# head -n 1 .DAJIN_temp/midsmask/barcode32_control.csv
# head -n 2 .DAJIN_temp/midsmask/tmp_maskDel
