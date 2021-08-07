# ちいさめ欠失
Read name = 0c061bae-a8b6-4a1e-80bb-a2330419283c
# おおきめ欠失
Read name = ae1eaf14-42cd-4e89-b5e3-5086b4d698dc
Read name = 62eb3fad-590c-4bb6-88a9-8c90c9515196
#
cat .DAJIN_temp/sam/"$sample_name"_control.sam | grep a2330419283c
cat .DAJIN_temp/sam/"$sample_name"_control.sam | grep 5086b4d698dc

cat .DAJIN_temp/mids/"$sample_name"_control.csv | grep a2330419283c | grep D
cat .DAJIN_temp/mids/"$sample_name"_control.csv | grep 5086b4d698dc | grep D

cat .DAJIN_temp/midsmask/"$sample_name"_control.csv | grep a2330419283c | grep D
cat .DAJIN_temp/midsmask/"$sample_name"_control.csv | grep 5086b4d698dc | grep D

cat .DAJIN_temp/score/barcode25_control.csv | grep -e a2330419283c -e 5086b4d698dc -e 8c90c9515196 |
  awk -F, '{for(i=2;i<=NF;i++) sum+=$i; print $1, sum}'

cat .DAJIN_temp/scalar/barcode25_control.csv | grep -e a2330419283c -e 5086b4d698dc

# 小さめと大きめの２アレルを１０本ずつ分けてみる
cat .DAJIN_temp/score/barcode25_control.csv | cut -d, -f1 | sort -u >tmp_id

cat .DAJIN_temp/score/barcode25_control.csv |
  awk -F, 'BEGIN {OFS=","} {
    sum=0
    for(i=(2+NF/3*1); i<=(2+NF/3*2); i++)
      sum+=$i

    if (sum < 4000000)
      $1="small-del"
    else
      $1="large-del"
    }1' >tmp_del

cat tmp_del | cut -d, -f1 >.DAJIN_temp/clustering/tmp_id.csv
cat tmp_del | cut -d, -f2- >.DAJIN_temp/clustering/tmp_score.csv
# grep small tmp_del | head -n 1000 | cut -d, -f1 >.DAJIN_temp/clustering/tmp_id.csv
# grep large tmp_del | head -n 1000 | cut -d, -f1 >>.DAJIN_temp/clustering/tmp_id.csv
# grep small tmp_del | head -n 1000 | cut -d, -f2- >.DAJIN_temp/clustering/tmp_score.csv
# grep large tmp_del | head -n 1000 | cut -d, -f2- >>.DAJIN_temp/clustering/tmp_score.csv

python .DAJIN_temp/library/clustering.py .DAJIN_temp/clustering/tmp_score.csv "$threads" |
  paste -d, - .DAJIN_temp/clustering/tmp_id.csv |
  cat >tmp

cat tmp |
  sed "s/^/${classif},/" |
  awk -F, 'BEGIN {OFS=","; clust_num=1} {
    allele=$1
    clust=$2
    id=$1","$2
    if(array_allele[allele]=="") {
      array_allele[allele]++
      clust_num=1
    }
    if(array_id[id]=="") {
      array_id[id]++
      allele_num[id]=clust_num
      clust_num++
      }
    $2=allele_num[id]
  }1' |
  sort | uniq -c

cat .DAJIN_temp/midsmask/barcode25_control.csv |
  awk -F, '{
    sum=0
    for(i=2; i<=NF; i++)
      if($i=="D")
        sum++
    print $1, sum}' |
  grep -e a2330419283c -e 5086b4d698dc -e 8c90c9515196
cat .DAJIN_temp/clustering/barcode25.csv | cut -d, -f1,2 | sort | uniq -c

cat .DAJIN_temp/midsmask/"$sample_name"_control.csv | grep 9995691e97b5 | grep D
cat .DAJIN_temp/midsmask/"$sample_name"_control.csv | grep 5086b4d698dc | grep D

#* メモ：STX2でcontrolとinversionにはないのにtargetのみある変なリード
# filterで取り除かれたよう
## 0689a48c-936a-47c9-b080-c2161d8f3683
