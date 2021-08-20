#!/bin/sh

maskDel="$(find .DAJIN_temp -name "maskDel.R")"

find .DAJIN_temp/midsmask/"$sample_name"* >.DAJIN_temp/midsmask/tmp_sample
find .DAJIN_temp/midsmask/"$control_name"* >.DAJIN_temp/midsmask/tmp_control

paste .DAJIN_temp/midsmask/tmp_sample .DAJIN_temp/midsmask/tmp_control |
  while read -r line; do
    set $line
    Rscript --vanilla --slave "$maskDel" "$1" "$2" \
      >.DAJIN_temp/midsmask/tmp_maskDel
    mv .DAJIN_temp/midsmask/tmp_maskDel "$1"
  done

# set .DAJIN_temp/midsmask/barcode25_control.csv .DAJIN_temp/midsmask/barcode30_control.csv
# cp $1 tmp_s
# cp $2 tmp_c
# set --
# cp .DAJIN_temp/midsmask/barcode32_control.csv tmp_c
# Rscript --vanilla --slave "$maskDel" tmp_s tmp_c >tmp_del
# head -n 1 .DAJIN_temp/midsmask/barcode31_control.csv
# head -n 1 .DAJIN_temp/midsmask/barcode32_control.csv
# head -n 2 .DAJIN_temp/midsmask/tmp_maskDel
