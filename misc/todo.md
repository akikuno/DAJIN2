# Todo

+ [ ] テストの用意
  + [ ] main
  + [ ] mids

## padding

+ paddingはmids-filterのときだけに使えば良い。
+ その他はmatchとして扱ってよさそ合う

## Ayabe-task1

+ `27db941cae33`はfloxなのにcontrolとされている
+ mids-maskのときに最大頻度の塩基でマスクする方法では、controlとfloxのヘテロのときにどっちか適当に割り当てられてしまうためにエラーとなる。
  + （案１）まずはfloxとcontrolをざっくりとわけて、次にマスク。そしてマスク塩基箇所はマスク塩基を除いたMIDS頻度に従ってマスク塩基をランダムにMIDSを付与すれば良いかも
    + controlにないloxpのなかにあるマスクは？
  + （案２）まずはfloxとcontrolをざっくりとわけて、次にマスク。knock-inだけの問題なので、knock-inアレルの場合はKI配列は無視する（マスクしない）
  + knock-inのときはマスクしない. それ以外はマスクする. 
    + → PMがヘテロのときには対応できない…
      + → マスクの仕方を、非マスク塩基のMIDS頻度にしたがってランダムに分配する方針にする
  + ~~（案３） マスク塩基はマスク塩基内のMIDS頻度にしたがって付与する？ <- マスク塩基内はクオリティが低いので却下~~

```bash
cat .DAJIN_temp/clustering/* | grep 27db941cae33
cat .DAJIN_temp/sam/barcode31_control.sam | grep 27db941cae33 | tr "[ACGT]" "," | grep [acgt]
cat .DAJIN_temp/sam/barcode31_flox.sam | grep 27db941cae33 | tr "[ACGT]" "," | grep [acgt]
cat .DAJIN_temp/mids/barcode31_control.csv | grep 27db941cae33 |tr -d "M" | grep -e [0-9][0-9]  -e D
cat .DAJIN_temp/mids/barcode31_flox.csv | grep 27db941cae33 | tr -d "M" | grep -e [0-9][0-9]  -e D
cat .DAJIN_temp/midsmask/barcode31_control.csv | grep 27db941cae33 | tr -d "M" | grep -e [0-9][0-9]  -e D
cat .DAJIN_temp/midsmask/barcode31_flox.csv | grep 27db941cae33 | tr -d "M" | grep -e [0-9][0-9] -e D

cat .DAJIN_temp/midsmask/barcode31_flox.csv | grep 27db941cae33 | tr -d "M" | awk '{print gsub("D","")}'

grep -e "^@" -e 27db941cae33 -e 1996ee7e907b .DAJIN_temp/sam/barcode31_flox.sam >tmp.sam
grep -e "^@" -e 27db941cae33 -e 1996ee7e907b .DAJIN_temp/midsfilter/barcode31_flox.csv >tmp.csv
set tmp.sam tmp.csv
cat "$1" |
  fmtPhred |
  join -t, - "$2" |
  grep 27db941cae33 >tmp
maskMIDS $1 $2 > tmp
cat tmp | grep 27db941cae33 | awk '{print gsub("D","")}'
cat tmp | awk -F, '{for(i=3;i<=NF;i++) if($i=="D") sumD++; print length($2), NF-2, sumD}'
cat tmp | maskMS | awk '{print gsub("D","")}'
cat tmp | maskMS | Rscript --vanilla "$maskMIDS" | cut -d,

cat .DAJIN_temp/scalar/barcode31_control.csv | grep 27db941cae33
cat .DAJIN_temp/scalar/barcode31_flox.csv | grep 27db941cae33
cat .DAJIN_temp/score/barcode31_control.csv | grep 27db941cae33 >tmp1.csv
cat .DAJIN_temp/score/barcode31_flox.csv | grep 27db941cae33 >tmp2.csv

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
