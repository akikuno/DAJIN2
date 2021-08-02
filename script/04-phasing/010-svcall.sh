#!/bin/sh

#----------------------------------------------------------
echo "$(date +'%Y-%m-%d %H:%M:%S') SV allele detetction" >>log_DAJIN.txt
#----------------------------------------------------------

mkdir -p .DAJIN_temp/sv

. .DAJIN_temp/library/calcHotelling.sh

cat .DAJIN_temp/scalar/"$control_name"* |
  cut -d, -f2 |
  calcHotelling |
  awk -F, '$1 < 3.841459' | #qchisq(0.95,1)
  cat >.DAJIN_temp/sv/tmp_control

# LOFだと値が小さすぎる（きれいにマッピングされすぎる）リードもSVと判定されてしまうため,
# 中央値よりも値が小さいリードはすべて正常と判定させる.
median_score=$(
  sort -n .DAJIN_temp/sv/tmp_control |
    awk -v wc="$(wc -l <.DAJIN_temp/sv/tmp_control)" 'NR==int(wc/2)'
)

find .DAJIN_temp/scalar/"$sample_name"* |
  while read -r line; do

  done

cat .DAJIN_temp/sv/"$sample_name"_id_score.csv |
  sort -t, |
  join -t "," - .DAJIN_temp/sv/allele_length |
  cut -d "," -f 3- |
  awk -F, '{sum=0; for(i=1; i<=NF-1; i++) sum+=$i; print log(sum/$NF)}' |
  cat >.DAJIN_temp/sv/tmp_sample

cat .DAJIN_temp/sv/tmp_sample |
  awk -v median="$median" '{
    ($0<median) ? $0="normal" : $0="SV_candidate"
  }1' >.DAJIN_temp/sv/tmp_sample_median

python .DAJIN_temp/library/lof_novelty.py -t "$threads" \
  -c .DAJIN_temp/sv/tmp_control -s .DAJIN_temp/sv/tmp_sample |
  awk '{($1 == 1) ? $1="normal" : $1="SV"}1' |
  paste - .DAJIN_temp/sv/tmp_sample_median |
  awk '{if($1=="SV" && $2=="normal") $1="normal"; print $1}' |
  paste -d "," - .DAJIN_temp/sv/"$sample_name"_id_score.csv |
  cut -d, -f 1-3 |
  sed "s/^/$sample_name,/" |
  awk -F, 'BEGIN{OFS=","} {print $2,$3,$4,$1}' |
  sort >.DAJIN_temp/sv/anomaly_id_"$sample_name".csv
