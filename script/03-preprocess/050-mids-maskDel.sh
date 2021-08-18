#!/bin/sh

maskDel="$(find .DAJIN_temp -name "maskDel.R")"

find .DAJIN_temp/midsmask/"$sample_name"* >.DAJIN_temp/midsmask/tmp_sample
find .DAJIN_temp/midsmask/"$control_name"* >.DAJIN_temp/midsmask/tmp_control

paste .DAJIN_temp/midsmask/tmp_sample .DAJIN_temp/midsmask/tmp_control |
  while read -r sample control; do
    Rscript --vanilla --slave "$maskDel" "$sample" "$control" \
      >.DAJIN_temp/midsmask/tmp_maskDel
    mv .DAJIN_temp/midsmask/tmp_maskDel "$sample"
  done

# cp .DAJIN_temp/midsmask/barcode31_control.csv tmp_s
# cp .DAJIN_temp/midsmask/barcode32_control.csv tmp_c
# Rscript --vanilla --slave "$maskDel" tmp_s tmp_c >tmp_del
# head -n 1 .DAJIN_temp/midsmask/barcode31_control.csv
# head -n 1 .DAJIN_temp/midsmask/barcode32_control.csv
# head -n 2 .DAJIN_temp/midsmask/tmp_maskDel
