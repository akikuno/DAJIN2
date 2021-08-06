#!/bin/sh

#----------------------------------------------------------
timestamp "SV allele detetction" >>log_DAJIN.txt
#----------------------------------------------------------

mkdir -p .DAJIN_temp/sv/
. .DAJIN_temp/library/calcHotelling.sh

cat .DAJIN_temp/classif/"$sample_name".csv |
  tee .DAJIN_temp/sv/tmp_sample |
  cut -d, -f3 |
  cat >.DAJIN_temp/sv/tmp_sample_score

control_scalar=$(find .DAJIN_temp/scalar/"$control_name"*)
# LOFだと値が小さすぎる（きれいにマッピングされすぎる）リードもSVと判定されてしまうため,
# 中央値よりも値が小さいリードはすべて正常と判定させる.
median_score=$(
  cat $control_scalar |
    awk -F, '$NF ~ /e-/ {$NF=0}1' |
    sort -t, -n -k2,2 |
    awk -F, -v wc="$(wc -l <$control_scalar)" 'NR==int(wc/2) {print $NF}'
)

cat "$control_scalar" |
  cut -d, -f2 |
  calcHotelling |
  paste -d, "$control_scalar" - |
  awk -F, -v median="$median_score" '!($2 > median && $NF > 3.841459) {print $2}' | #qchisq(0.95,1)
  cat >.DAJIN_temp/sv/tmp_control_score

python .DAJIN_temp/library/svLof.py \
  -c .DAJIN_temp/sv/tmp_control_score \
  -s .DAJIN_temp/sv/tmp_sample_score \
  -t "$threads" |
  awk '{($1 == 1) ? $1="normal" : $1="SV"}1' |
  paste -d, - .DAJIN_temp/sv/tmp_sample .DAJIN_temp/sv/tmp_sample_score |
  awk -F, -v median="$median_score" '{
    if($1=="SV" && $NF<median) $1="normal"
    if($1=="SV") $2="SV"
    print $2","$3
   }' |
  cat >.DAJIN_temp/sv/"$sample_name".csv

rm .DAJIN_temp/sv/tmp*

#? DEBUG ===
# head .DAJIN_temp/sv/"$sample_name".csv
# cat .DAJIN_temp/sv/"$sample_name".csv | cut -d, -f1 | sort | uniq -c
# cat .DAJIN_temp/sam/"$sample_name"_target.sam | grep 3284ac08-c086-4cbc-b871-8ddab98f42a0 | grep "[acgt]"
# cat .DAJIN_temp/sam/"$sample_name"_wt.sam | grep 3284ac08-c086-4cbc-b871-8ddab98f42a0 | grep "[acgt]"

# cat .DAJIN_temp/scalar/"$sample_name"_target.csv | grep 3284ac08-c086-4cbc-b871-8ddab98f42a0
# cat .DAJIN_temp/scalar/"$sample_name"_wt.csv | grep 3284ac08-c086-4cbc-b871-8ddab98f42a0
