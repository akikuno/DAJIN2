#!/bin/sh
# shellcheck disable=SC1091,SC2120

################################################################################
# input: MIDS (id, MIDS)
# output: MIDS score (id, I, D, S)
################################################################################

expansion() {
  if [ -p /dev/stdin ] && [ _"$*" = _"" ]; then
    cat -
  elif [ -r "$1" ]; then
    cat "$1"
  else
    echo "$*"
  fi |
    awk -F, '
    function mids_score(I,D,S) {
      $(i+(LEN)*0)=I
      $(i+(LEN)*1)=D
      $(i+(LEN)*2)=S
    }
    BEGIN {OFS=","} {
    LEN=NF-1
    for(i=2;i<=NF;i++) {
      if($i=="M")
        mids_score(0,0,0)
      else if ($i~/[0-9]+M/) {
        sub("M$","",$i)
        mids_score($i,0,0)
        }
      else if ($i~/[0-9]+D/) {
        sub("D$","",$i)
        mids_score($i,1,0)
      }
      else if ($i~/[0-9]+S/) {
        sub("S$","",$i)
        mids_score($i,0,1)
      }
      else if($i=="D") {
        mids_score(0,1,0)
      }
      else if($i=="S") {
        mids_score(0,0,1)
      }
    }
  }1'
}

rowSums() {
  if [ -p /dev/stdin ] && [ _"$*" = _"" ]; then
    cat -
  elif [ -r "$1" ]; then
    cat "$1"
  else
    echo "$*"
  fi |
    cat tmp.csv |
    awk -F, 'BEGIN {OFS=","} {
    for(i=2;i<=NF;i++) {
      if($i>0 && $(i+1)>0) {
        $i="@"
      }
      if($(i-1)=="@" && $i==1) {
        $i="@"
      }
    }
  }1' |
    sed "s/@,@/@@/g" |
    awk -F, 'BEGIN {OFS=","} {
    for(i=2; i<=NF; i++) {
      _rep=""
      if($i~/@/) {
        n=gsub("@",$i)
        for(j=1; j<=n; j++) {
          _rep=_rep n ","
        }
        gsub(",$" ,"" , _rep)
        $i=_rep
      }
    }
  }1'
}

midsToscore() {
  if [ -p /dev/stdin ] && [ _"$*" = _"" ]; then
    cat -
  elif [ -r "$1" ]; then
    cat "$1"
  else
    echo "$*"
  fi |
    expansion >.DAJIN_temp/tmp_expansion
  rowSums .DAJIN_temp/tmp_expansion >.DAJIN_temp/tmp_row
  colSums .DAJIN_temp/tmp_expansion >.DAJIN_temp/tmp_col
  rowColSums .DAJIN_temp/tmp_row .DAJIN_temp/tmp_col
}
