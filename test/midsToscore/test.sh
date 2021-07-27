#!/bin/sh

. library/midsToscore.sh

#echo "aaa,M,M,1M,M" |
echo "aaa,M,1S,D,M" |
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
      else if ($i~/[0-9]+[DSV]/) {
        sub("[DSV]$","",$i)
        mids_score($i+1,0,0)
      }
      else if($i=="D") {
        mids_score(0,1,0)
      }
      else if($i=="S") {
        mids_score(0,0,1)
      }
    }
  }1'

echo "aaa,0,0,1,0, 0,0,0,0, 0,0,0,0, 0,0,0,0"
