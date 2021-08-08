#!/bin/sh

################################################################################
# Input: CSV including allele, MIDS (e.g: normal_control,M,M,M,M,M,M,M)
# Output: CSV including allele, loc, M freq, I freq, D freq, and S freq (e.g: normal_control,998,0.888889,0,0.111111,0)
################################################################################

calcFreqMIDS() (
  # Input from pipe or file or or stdin
  if [ -p /dev/stdin ] && [ _"$*" = _"" ]; then
    cat -
  elif [ -r "$*" ]; then
    cat "$*"
  else
    echo "$*"
  fi |
    awk -F, 'BEGIN {OFS=","} {
    num_of[$1]++
    nf_of[$1]=NF
    for(i=2;i<=NF;i++) {
      mids[$1,i,$i]++
      }
    } END {
    for(type in num_of) {
      for(i=2;i<=nf_of[type];i++) {
        (mids[type, i, "M"]) ? M=mids[type, i, "M"]/num_of[type] : M=0
        (mids[type, i, "I"]) ? I=mids[type, i, "I"]/num_of[type] : I=0
        (mids[type, i, "D"]) ? D=mids[type, i, "D"]/num_of[type] : D=0
        (mids[type, i, "S"]) ? S=mids[type, i, "S"]/num_of[type] : S=0
        print type,i-1,M,I,D,S
      }
    }
  }'
)
