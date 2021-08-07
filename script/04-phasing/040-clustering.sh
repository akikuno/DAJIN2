#!/bin/sh

timestamp "Clustering" | tee log_DAJIN.txt

mkdir -p .DAJIN_temp/clustering/

cat .DAJIN_temp/sv/"$sample_name".csv |
  cut -d, -f 1 |
  sort |
  uniq -c |
  while read -r num classif; do
    allele="$classif"
    [ _"$classif" = "_SV" ] && allele=control
    #
    grep "^$classif" .DAJIN_temp/sv/"$sample_name".csv |
      cut -d, -f2 |
      sort -u >.DAJIN_temp/clustering/tmp_id.csv
    #
    if [ "$num" -gt 5 ]; then
      cat .DAJIN_temp/score/"$sample_name"_"$allele".csv |
        join -t, - .DAJIN_temp/clustering/tmp_id.csv |
        cut -d, -f2- >.DAJIN_temp/clustering/tmp_score.csv
      python .DAJIN_temp/library/clustering.py .DAJIN_temp/clustering/tmp_score.csv "$threads"
    else
      awk -v num="$num" 'BEGIN{for(i=1;i<=num;i++) print 1}'
    fi |
      paste -d, - .DAJIN_temp/clustering/tmp_id.csv |
      sed "s/^/${classif},/"
  done |
  #* format allele numbers
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
  sort -u >.DAJIN_temp/clustering/"$sample_name".csv

rm .DAJIN_temp/clustering/tmp* 2>/dev/null || :
rm -rf .DAJIN_temp/clustering/joblib/ 2>/dev/null || :
# ###? DEBUG
# cat .DAJIN_temp/clustering/"$sample_name".csv |
#   cut -d, -f 1-2 |
#   sort |
#   uniq -c

# echo "CLUSTERING FINISH"
# exit 0
