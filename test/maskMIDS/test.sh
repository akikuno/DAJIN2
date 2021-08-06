#!/bin/sh

. library/samToMIDS.sh

minimap2 -ax map-ont --cs=long test/maskMIDS/ref.fa test/maskMIDS/que.fq >tmp.sam
minimap2 -ax map-ont --cs=long test/maskMIDS/ref_long.fa test/maskMIDS/que_long.fq >tmp_long.sam

cat tmp.sam | samToMIDS >tmp_mids.csv
cat tmp_long.sam | samToMIDS >tmp_mids_long.csv
cat tmp_mids.csv

SAM=tmp.sam
MIDS=tmp_mids.csv

SAM=tmp_long.sam
MIDS=tmp_mids_long.csv

set "$SAM" "$MIDS"
set .DAJIN_temp/sam/barcode31_albino.sam .DAJIN_temp/mids/barcode31_albino.csv
set .DAJIN_temp/sam/barcode31_control.sam .DAJIN_temp/mids/barcode31_control.csv
maskMIDS "$SAM" "$MIDS"

cat "$SAM" |
  grep -v "^@" |
  sort -k 1,1 -k 4,4n |
  cut -f 1,6,11 |
  awk '
  function trimcrip(cigar, score) {
    match(cigar, /[0-9][0-9]*S/)
    leftclip=substr(cigar, RSTART, RLENGTH-1)
    sub(/[0-9][0-9]*S/, "", cigar)
    match(cigar, /[0-9][0-9]*S/)
    rightclip=substr(cigar, RSTART, RLENGTH-1)
    score=substr(score, leftclip+1)
    score=substr(score, 1, length(score)-rightclip)
    return score
  }
  {
    n_id[$1]++
    cigar=$2
    score=$3
  if(cigar~/S/) {
    score_of[$1] = score_of[$1] trimcrip(cigar, score)
  }
  else {
    score_of[$1] = score_of[$1] score
  }
  } END {
    for(id in score_of) {
      print id","score_of[id]
    }
  }' |
  join -t, - "$MIDS" |
  awk -F, 'BEGIN {OFS=","} {
    readid=$1
    split($2, cigar, "")
    idx=1
    for (i=3; i<=NF; i++) {
      if ($i ~ /[MS]/ && cigar[idx] ~ /[!"#$%&'\''()*]/) {
        $i="N"
        idx++
      }
      else if ($i == "D") {
        idx=idx
      }
      else if ($i ~ /[0-9]/) {
        I=$i
        sub(/[MDSmis]$/, "", I)
        idx=idx+I+1
      }
      else {
        idx++
      }
    }
  }1' |
  cut -d, -f1,3- |
  ./test_table.R

# maskMIDS
## 最大頻度のMIDSを抽出
## Nを最大頻度のMIDSに置換する
## すべてNの場合はMに置換する

cp tmp_ tmp.csv
cat tmp.csv | ./test_table.R

cat tmp_ |
  cut -d, -f1-3 |
  awk -F, '{
    id[NR]=$1
    for (i=2; i<=NF; i++) {
      mids[i, NR]=$i
      count[i,$i]++
    }
  } END {
    for (key in count) {
      split(key, sep, SUBSEP)
      num=count[key]
      i=sep[1]
      midsn=sep[2]
      if (num==NR && midsn=="N") {
        for (nr=1;nr<=NR;nr++){
          mids[i, nr] = "M"
        }
      }
      else if (num==NR) {
        for (nr=1;nr<=NR;nr++){
          mids[i, nr] = midsn
        }
      }
      else if (num<NR) {
        for (nr=1;nr<=NR;nr++){
          mids[i, nr] = midsn
        }
      }
    }
    # output
    for (nr=1; nr<=NR; nr++) {
      printf id[nr]","
      for(i=2; i<=NF; i++) {
        printf mids[i, nr]","
      }
      print ""
    }
  }'
