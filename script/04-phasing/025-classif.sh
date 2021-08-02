#!/bin/sh

#----------------------------------------------------------
timestamp "Allele classification" >>log_DAJIN.txt
#----------------------------------------------------------

mkdir -p .DAJIN_temp/classif

find .DAJIN_temp/scalar/"$sample_name"* |
  while read -r line; do
    allele=$(basename "${line%.csv}" | cut -d_ -f2-)
    awk -v al="$allele" '{print al","$0}' "$line"
  done |
  awk -F, '{
    allele=$1
    read=$2
    score=$3
    if(score_of[read] == "") score_of[read]="inf"
    if(score_of[read]>score) {
      output[read]=$0
      score_of[read]=score
    }} END {
      for(read in output) print output[read]
    }' |
  cat >.DAJIN_temp/classif/"$sample_name".csv
