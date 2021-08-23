#!/bin/sh

#----------------------------------------------------------
timestamp "SV allele detetction" log_DAJIN.txt
#----------------------------------------------------------

mkdir -p .DAJIN_temp/sv/

# control =================================================

# control_scalarのSVがあるとsampleのSV判定が甘くなるので、事前にcontrol_scalarのSVを取り除きます。
cat .DAJIN_temp/midsmask/"$control_name"*_control.csv |
  sed "s/,/ /" |
  awk '{
    numI = gsub(/[5-9][0-9][0-9]*/, "@", $2)
    numD = gsub("(D,){50}", "@", $2)
    numS = gsub("(S,){50}", "@", $2)
    if (numI + numD + numS == 0)
      print $1
  }' |
  sort |
  join -t, - .DAJIN_temp/scalar/"$control_name"_control.csv |
  cat >.DAJIN_temp/sv/tmp_control_scalar

# LOFだと値が小さすぎる（きれいにマッピングされすぎる）リードもSVと判定されてしまうため,
# 中央値よりも値が小さいリードはすべて正常と判定させる.
median_score=$(
  cat .DAJIN_temp/sv/tmp_control_scalar |
    awk -F, '$NF ~ /e-/ {$NF=0}1' |
    sort -t, -n -k2,2 |
    awk -F, -v wc="$(wc -l <.DAJIN_temp/sv/tmp_control_scalar)" 'NR==int(wc/2) {print $NF}'
)

cat .DAJIN_temp/sv/tmp_control_scalar |
  cut -d, -f2 |
  calcHotelling |
  paste -d, .DAJIN_temp/sv/tmp_control_scalar - |
  awk -F, -v median="$median_score" '!($2 > median && $NF > 3.841459) {print $2}' | #qchisq(0.95,1)
  cat >.DAJIN_temp/sv/tmp_control_score

# sample =================================================

cat .DAJIN_temp/classif/"$sample_name".csv |
  tee .DAJIN_temp/sv/tmp_sample |
  cut -d, -f3 |
  cat >.DAJIN_temp/sv/tmp_sample_score

# SV detection  =================================================
svLof="$(find .DAJIN_temp -name "svLof.py")"
python "$svLof" \
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
