################################################################################
cat <<EOF >>log_DAJIN.txt
==========================================================
Classify alleles...
==========================================================
EOF
################################################################################

mkdir -p .DAJIN_temp/classif

find .DAJIN_temp/score/"$sample_name"_* |
  while read -r line; do
    cat <<____EOF
    awk -F, -v sample="${sample_name}" 'BEGIN {OFS=","} {
      sum=0
      allele=FILENAME
      sub(".*"sample"_", "", allele)
      sub(".csv$", "", allele)
      for(i=2;i<=NF;i++) sum+=\$i
      print sum, allele, \$0
    }' "$line" > "${line%.csv}"_tmp &
____EOF
  done |
  awk '1; END {print "wait"}' |
  sh

cat .DAJIN_temp/score/"$sample_name"*_tmp |
  awk -F, '{
    score=$1; allele=$2; id=$3
    if(score_of[id] == "") score_of[id]="inf"
    if(score_of[id]>score) {
      score_of[id]=score; seq[id]=$0
    }} END {for(id in seq) print seq[id]}' |
  cut -d "," -f 2- |
  sort -t, |
  cat >.DAJIN_temp/classif/"$sample_name"_id_score.csv

rm .DAJIN_temp/score/"$sample_name"*_tmp
