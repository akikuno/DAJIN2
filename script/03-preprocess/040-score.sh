#!/bin/sh

#----------------------------------------------------------
echo "$(date +'%Y-%m-%d %H:%M:%S') MIDSV scoring" >>log_DAJIN.txt
#----------------------------------------------------------

mkdir -p .DAJIN_temp/score

head .DAJIN_temp/midsv/barcode31_albino.csv |
  grep 00328905-1c46-4f17-8816-7881d8d44bb3 |
  sed -e "s/,\=/,D/g" |
  sed -e "s/,M/,0/g" -e s"/,M$/,0/" |
  cmd='. .DAJIN_temp/library/samTomidsv.sh; samTomidsv '
