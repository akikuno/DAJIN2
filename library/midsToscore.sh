#!/bin/sh
# shellcheck disable=SC1091,SC2120

################################################################################
# input: MIDS (id, MIDS)
# output: MIDS score (id, I, D, S)
################################################################################

expansion() {
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

mutToat() {
  awk '
    BEGIN {
      FS=","
      OFS=""
    } {
    $1=$1","
    for (i=2; i<=NF; i++) {
      if ($i>0 && $(i+1)>0)
        $i="@"
      else if ($(i-1)=="@" && $i==1)
        $i="@,"
      else
        $i=$i","
    }
    sub(",$", "", $NF)
  }1'
}

atToscore() {
  awk '
    BEGIN {
      FS=","
      OFS=","
    } {
    for (i=2; i<=NF; i++) {
      n=""
      if ($i ~ /@/) {
        num=length($i)
        for (j=1; j<=num; j++) n=n num ","
        sub(",$", "", n)
        $i=n
      }
    }
  }1'
}

rowScore() {
  cat "$1" |
    mutToat |
    atToscore
}

colScore() {
  Rscript --vanilla .DAJIN_temp/library/colTable.R "$1"
}

rowColSums() {
  Rscript --vanilla .DAJIN_temp/library/rowColSums.R "$1" "$2"
}

midsToscore() {
  mkdir -p .DAJIN_temp
  if [ -p /dev/stdin ] && [ _"$*" = _"" ]; then
    cat -
  elif [ -r "$1" ]; then
    cat "$1"
  else
    echo "$*"
  fi |
    expansion >.DAJIN_temp/tmp_expansion
  rowScore .DAJIN_temp/tmp_expansion >.DAJIN_temp/tmp_row
  colScore .DAJIN_temp/tmp_expansion >.DAJIN_temp/tmp_col
  rowColSums .DAJIN_temp/tmp_row .DAJIN_temp/tmp_col
  rm .DAJIN_temp/tmp_*
}
