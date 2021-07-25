#!/bin/sh

. library/midsv_scoring.sh

echo "aaa,M,M,1M,M" |
  awk -F, '
  function midsv_score(I,D,S,V) {
    $(i+(LEN-1)*0)=I
    $(i+(LEN-1)*1)=D
    $(i+(LEN-1)*2)=S
    $(i+(LEN-1)*3)=V
  } BEGIN {OFS=","} {
    LEN=NF-1
    for(i=2;i<=NF;i++) {
      if($i=="M")
        midsv_score(0,0,0,0)
      else if ($i~/[0-9]+M/) {
        sub("M$","",$i)
        midsv_score($i,0,0,0)
        }
      else if ($i~/[0-9]+[DSV]/) {
        sub("[DSV]$","",$i)
        midsv_score($i+1,0,0,0)
      }
      else if($i=="D")
        midsv_score(0,1,0,0)
      else if($i=="S")
        midsv_score(0,0,1,0)
      else if($i=="V")
        midsv_score(0,0,0,1)
    }
  }1'

echo "aaa,0,0,1,0, 0,0,0,0, 0,0,0,0, 0,0,0,0"
