#!/bin/sh

fmtScore() {
  cat - |
    grep -v "^@" |
    grep "cs:Z" |
    sort -k 1,1 -k 4,4n |
    cut -f 1,6,11 |
    awk '{gsub(",", "_", $NF)}1' |
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
      } {
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
    sort -t,
}

maskMS() {
  awk -F, 'BEGIN {OFS=","} {
    readid=$1
    split($2, cigar, "")
    idx=1
    for (i=3; i<=NF; i++) {
      if ($i ~ /[0-9]/) {
        I=$i
        sub(/[MDSmds]$/, "", I)
        idx=idx+I+1
      }
      else if ($i == "D") {
        idx=idx
      }
      else if ($i ~ /[MS]/ && cigar[idx] ~ /[!"#$%&'\''()*]/) {
        $i="N"
        idx++
      }
      else {
        idx++
      }
    }
  }1' |
    cut -d, -f1,3-
}

maskMIDS() {
  maskMIDS="$(find .DAJIN_temp/ -name "maskMIDS.R")"
  cat "$1" |
    fmtScore |
    join -t, - "$2" |
    maskMS |
    Rscript --vanilla "$maskMIDS"
}
