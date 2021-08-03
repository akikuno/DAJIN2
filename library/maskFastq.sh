#!/bin/sh

maskFastq() {
  if [ -p /dev/stdin ] && [ _"$*" = _"" ]; then
    cat -
  elif [ -r "$*" ]; then
    cat "$*"
  else
    echo "$*"
  fi |
    sed "s/^@/\n@/" |
    awk 'BEGIN{RS=""; FS="\n"; OFS="\n"} {
    split($2, sequence, "")
    n=split($NF, quality, "")
    for(i=1 && seq=""; i<=n; i++) {
      if (quality[i] ~ "[!\"#$%&'\''()*]") {
        seq=seq "N"
      }
      else {
        seq=seq sequence[i]
      }
    }
    $2=seq
  }1'
}
