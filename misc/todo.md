# Todo

## Ayabe-task1

+ `27db941cae33`はfloxなのにcontrolとされている
+ `205a6973c353`はfloxなのに片方のloxPがmaskされている midsmaskの問題
```bash
cat .DAJIN_temp/clustering/* | grep 205a6973c353
cat .DAJIN_temp/mids/barcode31_control.csv | grep 205a6973c353 | grep -e 45 -e 68
cat .DAJIN_temp/midsmask/barcode31_control.csv | grep 205a6973c353 | grep -e 45 -e 68

cat .DAJIN_temp/sam/barcode31_control.sam | grep -e "^@" -e 205a6973c353 >tmp.sam
cat .DAJIN_temp/mids/barcode31_control.csv | grep -e "^@" -e 205a6973c353 >tmp.csv
set tmp.sam tmp.csv
maskMIDS="$(find .DAJIN_temp/ -name "maskMIDS.R")"
cat "$1" |
  fmtScore |
  join -t, - "$2" |
  maskMS |
  Rscript --vanilla "$maskMIDS"
```

+ Cables2 flox knockinが floxではなくcontrolになっている.
  + `67803dec0a5e`はfloxなのにcontrolのほうがscoreが高くなっている.

```bash
cat .DAJIN_temp/scalar/barcode31_control.csv | grep 67803dec0a5e
cat .DAJIN_temp/scalar/barcode31_flox.csv | grep 67803dec0a5e
cat .DAJIN_temp/score/barcode31_control.csv | grep 67803dec0a5e | awk -F, '{while(i!=(NF-1)/3) print i,$(i++)}' > tmp
cat .DAJIN_temp/midsmask/barcode31_control.csv | grep 67803dec0a5e | grep [0-9][0-9] |
awk -F, '{for(i=2;i<=NF;i++) {if ($i~/[0-9][0-9]M/) {print i,$i}}}'
cat tmp | grep -e 1735 -e 2383

echo $((38*38)) $((62*62))
cat << EOF>tmp
hoge,5M,D,D,D,D,D
fuga,5M,D,D,D,D,D
EOF

echo "hoge,M,M,5M,M,62M,D,D,D,D,D" >tmp
echo "fuga,M,M,5M,M,M,D,D,D,D,D" >>tmp
. $(find library -name "midsToScore.sh")
midsToScore tmp
set tmp
expansion "$1" >tmp_expansion
rowScore tmp_expansion >tmp_row
colScore tmp_expansion >tmp_col
Rscript .DAJIN_temp/library/04-phasing/colScore.R tmp_expansion > tmp_col
cat tmp_row tmp_col
rowColSums tmp_row tmp_col
rowColMul tmp_row tmp_col
Rscript .DAJIN_temp/library/04-phasing/rowColMul.R tmp_row tmp_col

cat .DAJIN_temp/scalar/barcode31_flox.csv | grep 67803dec0a5e
cat .DAJIN_temp/scalar/barcode31_control.csv | grep 67803dec0a5e
```

## Conseusns
+ HTML/FASTAレポート

## アレル頻度のプロット
+ DAJIN2_tmp01の過去の遺産をコピー

## Clustering:ファジーな変異情報を無視するためにはどうしましょう？
+ PhredスコアとMIDS頻度情報を用いる
  + MIDS頻度を用いて
  + Phredスコアが低いものを, 高いものに合わせる
    + → albino部位の欠失で

# Done

## 2021-08-11
+ [x] ayabe-task1のfloxのサンプルを例題に追加

## 2021-07-26
+ [x] paddingの”=”のあつかいはどうしよう. MatchかDeletionか？ → ひとまずDeletionとしてあつかう
+ [x] inversionアレルの変異情報が消えてしまう → 小文字にすることで完成しました. 

---

### メモ
