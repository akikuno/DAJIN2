#!/bin/sh
# shellcheck disable=SC1091,SC2120

################################################################################
# input: MIDS (id, MIDS)
# output: MIDS score (id, I, D, S)
################################################################################

input() {
  if [ -p /dev/stdin ] && [ _"$*" = _"" ]; then
    cat -
  elif [ -r "$1" ]; then
    cat "$1"
  else
    echo "$*"
  fi
}

midsToscore() {
  input |
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
