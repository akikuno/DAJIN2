#!/bin/sh

#----------------------------------------------------------
timestamp "SV allele detetction" >>log_DAJIN.txt
#----------------------------------------------------------

mkdir -p .DAJIN_temp/sv/
. .DAJIN_temp/library/calcHotelling.sh

cat .DAJIN_temp/classif/"$sample_name".csv |
  tee .DAJIN_temp/sv/tmp_sample |
  cut -d, -f3 |
  calcHotelling >.DAJIN_temp/sv/tmp_sample_score

cat .DAJIN_temp/scalar/"$control_name"* |
  cut -d, -f2 |
  calcHotelling |
  awk -F, '$1 < 3.841459' | #qchisq(0.95,1)
  cat >.DAJIN_temp/sv/tmp_control_score

# LOFだと値が小さすぎる（きれいにマッピングされすぎる）リードもSVと判定されてしまうため,
# 中央値よりも値が小さいリードはすべて正常と判定させる.
median_score=$(
  sort -n .DAJIN_temp/sv/tmp_control_score |
    awk -v wc="$(wc -l <.DAJIN_temp/sv/tmp_control_score)" 'NR==int(wc/2)'
)

python .DAJIN_temp/library/svLof.py \
  -c .DAJIN_temp/sv/tmp_control_score \
  -s .DAJIN_temp/sv/tmp_sample_score \
  -t "$threads" |
  awk '{($1 == 1) ? $1="normal" : $1="SV"}1' |
  paste -d, - .DAJIN_temp/sv/tmp_sample |
  awk -F, -v median="$median_score" '{
    if($1=="SV" && $NF<median) $1="normal"
    if($1=="SV") $2="SV"
    print $2","$3
   }' |
  cat >.DAJIN_temp/sv/"$sample_name".csv

rm .DAJIN_temp/sv/tmp*

#? DEBUG ==========
# cat .DAJIN_temp/sv/"$sample_name".csv | cut -d, -f1 | sort | uniq -c
# cat .DAJIN_temp/sv/"$sample_name".csv | grep control, | cut -d, -f 2 | sort -u >tmpcontrol
# cat .DAJIN_temp/sv/"$sample_name".csv | grep albino, | cut -d, -f 2 | sort -u >tmpalbino

# cat .DAJIN_temp/classif/"$sample_name".csv | sort -t, -k2,2 | join -1 2 -2 1 -t, - tmpcontrol
# cat .DAJIN_temp/classif/"$sample_name".csv | sort -t, -k2,2 | join -1 2 -2 1 -t, - tmpalbino

# cat .DAJIN_temp/midsmask/barcode31_control.csv | sort -t, | join -t, - tmpcontrol |
#   awk -F, '{loc=829; print $(loc+1)}' | sort | uniq -c
