#!/bin/sh

#==============================================================================
# Read name = e5285832-50a3-4ef8-bf62-6150dc025be3がalbinoなのにcontrolとなっている
### 点変異の位置は829番目
### albinoポジコンは3533774adc91
#==============================================================================

que=6150dc025be3
ref=3533774adc91

cat .DAJIN_temp/clustering/*.csv | grep "$que"
cat .DAJIN_temp/clustering/*.csv | grep "$ref"
cat .DAJIN_temp/mids/"$sample_name"*_control.csv | grep "$que" | awk -F, '{print $(829+1)}'

grep "$que" .DAJIN_temp/scalar/*.csv
grep "$ref" .DAJIN_temp/scalar/*.csv

#==============================================================================
# きちんとSのスコアが反映されているのか確認する
#==============================================================================

grep "$ref" .DAJIN_temp/score/*.csv |
  # awk -F, '{for(i=1;i<=NF;i++) if($i==3753) print i}'
  awk -F, '{nf=NF+1; print $1, $(nf/3*2+829)}'

grep "$que" .DAJIN_temp/score/*.csv |
  awk -F, '{nf=NF+1; print $1, $(nf/3*2+829)}'

# 反映されている！！

#==============================================================================
#
#==============================================================================
grep "$ref" .DAJIN_temp/score/*.csv |
  awk -F, '{nf=NF+1; print $1, $(nf/3*2+829)}'

grep "$que" .DAJIN_temp/midsmask/*_control.csv
grep "$que" .DAJIN_temp/score/*_control.csv |
  awk -F, '{for(i=2;i<=NF;i++) if ($i>500) print i,$i}' | sort >tmp1
grep "$que" .DAJIN_temp/score/*_albino.csv |
  awk -F, '{for(i=2;i<=NF;i++) if ($i>500) print i,$i}' | sort >tmp2

join -v 1 tmp2 tmp1
idx=1000
idx=830
cat .DAJIN_temp/midsmaskByPhred/barcode31_control.csv |
  awk -F, -v i="$idx" '{print i,$i}' |
  sort | uniq -c
cat .DAJIN_temp/midsmaskByPhred/barcode32_control.csv |
  awk -F, -v i="$idx" '{print i,$i}' |
  sort | uniq -c

# 自然な（こじつけではない）補正法は？
# ControlからQueryを引いて、プラスマイナスが強いほう（絶対値が大きい）が興味のある変異である！
# 差の絶対値が0.1以下のものは、シークエンスエラーとみなして良さそう
# 差の絶対値が0.1以上のものは、真の変異の可能性が高い

# 988
grep "$que" .DAJIN_temp/score/*_control.csv |
  awk -F, '{i=988; print i,$i}'
grep "$que" .DAJIN_temp/score/*_albino.csv |
  awk -F, '{i=988; print i,$i}'

cat .DAJIN_temp/midsmaskByPhred/*_control.csv |
  awk -F, '{i=988; print i,$i}' |
  sort | uniq -c

cat .DAJIN_temp/midsmaskByPhred/*_albino.csv |
  awk -F, '{i=988; print i,$i}' |
  sort | uniq -c

cat tmp_query.csv |
  awk -F, '{i=830; print i,$i}' |
  sort | uniq -c
cat tmp_control.csv |
  awk -F, '{i=830; print i,$i}' |
  sort | uniq -c
