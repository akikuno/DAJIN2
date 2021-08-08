#!/bin/sh

echo TEST #?????????????????
exit 0    #?????????????????
. .DAJIN_temp/library/calcHotelling.sh
mkdir -p .DAJIN_temp/consensus /tmp/DAJIN/consensus
control_scalar=.DAJIN_temp/scalar/"$control_name"_control.csv

# control =================================================

cat .DAJIN_temp/midsmask/barcode32_control.csv |

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
  awk -F, -v median="$median_score" '!($2 > median && $NF > 3.841459) {print $1}' | #qchisq(0.95,1)
  sort -u |
  join -t, .DAJIN_temp/score/"$control_name"_control.csv - |
  cat >.DAJIN_temp/consensus/tmp_control_score
