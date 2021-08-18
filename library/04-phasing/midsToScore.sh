#!/bin/sh
# shellcheck disable=SC1091,SC2120

################################################################################
# input: MIDS (id, MIDS)
# output: MIDS score (id, I, D, S)
################################################################################

expansion() {
  if [ -p /dev/stdin ] && [ "$#" -eq 0 ]; then
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
    function invToD(data) {
      id=$1
      gsub(/[mids]/, "D", $0)
      $1=id
    }
    BEGIN {OFS=","} {
    LEN=NF-1
    invToD($0)
    for(i=2;i<=NF;i++) {
      if($i=="M")
        mids_score(0,0,0)
      #* insertion->match
      else if ($i~/[0-9]+M/) {
        sub("M$","",$i)
        $i = $i * $i
        mids_score($i,0,0)
        }
      #* insertion->deletion
      else if ($i~/[0-9]+D/) {
        sub("D$","",$i)
        $i = $i * $i
        mids_score($i,1,0)
      }
      #* insertion->substitution
      else if ($i~/[0-9]+S/) {
        sub("S$","",$i)
        $i = $i * $i
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

mutToAt() {
  awk '
    BEGIN {
      FS=","
      OFS=""
    } {
    $1=$1","
    ins_nf=(NF-1)/3
    for (i=2; i<=ins_nf; i++) {
      $i=$i","
    }
    for (i=ins_nf+1; i<=NF; i++) {
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

atToScore() {
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
    mutToAt |
    atToScore
}

colScore() {
  colScore="$(find .DAJIN_temp/ -name "colScore.R")"
  Rscript --vanilla "$colScore" "$1"
}

# rowColSums() {
#   rowColSums="$(find .DAJIN_temp/ -name "rowColSums.R")"
#   Rscript --vanilla "$rowColSums" "$1" "$2"
# }

rowColMul() {
  rowColMul="$(find .DAJIN_temp/ -name "rowColMul.R")"
  Rscript --vanilla "$rowColMul" "$1" "$2"
}

midsToScore() {
  mkdir -p .DAJIN_temp
  suffix="${1##*/}"
  expansion "$1" >.DAJIN_temp/tmp_expansion_"$suffix"
  rowScore .DAJIN_temp/tmp_expansion_"$suffix" >.DAJIN_temp/tmp_row_"$suffix"
  colScore .DAJIN_temp/tmp_expansion_"$suffix" >.DAJIN_temp/tmp_col_"$suffix"
  rowColMul .DAJIN_temp/tmp_row_"$suffix" .DAJIN_temp/tmp_col_"$suffix"
  rm .DAJIN_temp/tmp_*"$suffix"
}
