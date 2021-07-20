#!/bin/sh

calc_hotelling() (
  if [ -p /dev/stdin ] && [ _"$*" = _"" ]; then
    cat -
  elif [ -r "$*" ]; then
    cat "$*"
  else
    echo "$*"
  fi |
    awk '{
      data[NR]=$0
      sum+=$1
    } END {
      mean=sum/NR
      for(i=1;i<=NR;i++) {
        residual+=(data[i]-mean)^2
      }
      var=residual/NR
      for(i=1;i<=NR;i++) {
        hotelling=(data[i]-mean)^2/var
        print hotelling
      }
    }'
)
