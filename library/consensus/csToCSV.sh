#!/bin/sh
# shellcheck disable=SC1091,SC2120

################################################################################
# input: SAM file using `minimap2 -ax splice --cs=long`
# output: DNA sequence with MIDS separation (location, MIDS, Nucreotide)
################################################################################

fmtSam() {
  cat - |
    awk '
    /^@SQ/ {
      for(i=1;i<=NF;i++) {
        if($i ~ /^SN:/) allele=$i
        if($i ~ /^LN:/) reflen=$i
      }
      sub(/.*SN:/,"", allele)
      sub(/.*LN:/,"", reflen)
    }
    /cs:Z:=/ {
      id=$1; flag=$2; start=$4
      for(i=1;i<=NF;i++) if ($i ~ /^cs:Z:=/) cstag=$i
      sub("cs:Z:=","",cstag)
      gsub("=", " ", cstag)
      gsub(/\+/, " +", cstag)
      gsub(/\-/, " -", cstag)
      gsub(/\*/, " *", cstag)
      gsub("~", " ~", cstag)
      $1=id","flag","start","allele","reflen
      $2=cstag
      print $1,$2
    }' |
    sort -t "," -k 1,1 -k 3,3n
}

fmtIDS() {
  awk '{
    for (i=2; i<=NF; i++) {
      if ($i ~ /^\-/) {
        $i = substr($i,2)
        n = split($i, a, "")
        seq = ""
        for (j=1; j<=n; j++) {
          seq = seq "D"
        }
        $i = seq
      }
      else if ($i ~ /^\*/) {
        sub(/\*[acgt]/, "", $i)
      }
      else if ($i ~ /^\+/) {
        $i = substr($i,2)
        $i = $i substr($(i+1),1,1)
        $(i+1) = substr($(i+1),2)
      }
    }
  }1'
}

padding() {
  awk '{
    start=$2-1
    end=$3
    len=NF-3+start
    pad_start=""
    pad_end=""
    for(i=1;i<=start;i++) pad_start=pad_start "D "
    $4=pad_start $4
    for(i=len;i<end;i++) pad_end=pad_end " D"
    $NF=$NF pad_end
  }1' |
    sed "s/  */ /g"
}

spaceTocomma() {
  sed -e "s/  */,/g" -e "s/,$//"
}

csToCSV() (
  if [ -p /dev/stdin ] && [ _"$*" = _"" ]; then
    cat -
  elif [ -r "$*" ]; then
    cat "$*"
  else
    echo "$*"
  fi |
    fmtSam |
    padding |
    spaceTocomma |
    awk -F, 'NF==$3+3' |
    cut -d, -f 1,4- |
    grep -v "^$" |
    sort -t,
)
