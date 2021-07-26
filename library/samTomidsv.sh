#!/bin/sh
# shellcheck disable=SC1091,SC2120

################################################################################
# input: SAM file using `minimap2 -ax splice --cs=long`
# output: Match(M), Insertion(1M), Deletion(D), Substitution(S), Inversion (V) and "= (null)"
################################################################################

fmt_sam() {
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
      gsub(/\+/, " ", cstag)
      gsub(/\-/, " -", cstag)
      gsub(/\*/, " *", cstag)
      gsub("~", " ~", cstag)
      $1=id","flag","start","allele","reflen
      $2=cstag
      print $1,$2
    }' |
    sort -t "," -k 3,3n
}

matchToM() {
  awk '{for(i=2;i<=NF;i++) gsub(/[ACGT]/, "M ", $i)}1'
}

subToS() {
  awk '{for(i=2;i<=NF;i++) gsub(/\*[acgt][acgt]/, "S ", $i)}1'
}

delToD() {
  awk '{
      for(i=2;i<=NF;i++) {
        if($i ~ /^-/) {
          str=""
          for(j=1;j<=int(length($i)-1);j++) str=str "D "
          $i=str
        }
      }}1'
}

insToNum() {
  awk '{
    for(i=3;i<=NF;i++) {
      if($i~/^[acgt]/){
        $i=length($i) $(i+1)
        $(i+1)=""
      }
    }}1'
}

padding() {
  awk '{
    start=$2-1
    end=$3
    len=NF-3+start
    pad_start=""
    pad_end=""
    for(i=1;i<=start;i++) pad_start=pad_start "= "
    $4=pad_start $4
    for(i=len;i<end;i++) pad_end=pad_end " ="
    $NF=$NF pad_end
  }1' |
    sed "s/  */ /g"
}

spaceTocomma() {
  sed -e "s/  */,/g" -e "s/,$//"
}

samTomidsv() (
  if [ -p /dev/stdin ] && [ _"$*" = _"" ]; then
    cat -
  elif [ -r "$1" ]; then
    cat "$1"
  else
    echo "$*"
  fi |
    # cat test/samTomidsv/input-del_to_wt_allele.sam |
    # cat test/samTomidsv/test_inv.sam |
    # cat test/samTomidsv/test_del_long.sam |
    fmt_sam |
    matchToM |
    subToS |
    delToD |
    awk '{$1=$1","}1' |
    #* Large deletion and Inversion -------------------------
    awk -F, '
      function padD(iter,    i,str) {
        for (i=1; i<=iter; i++) str=str " D "
        return str
      }

      function rm_insertion(string) {
        gsub("[acgt][acgt]*", "", string)
        return string
      }

      function csCat(c_of, s_of, iter,    i,cs) {
        for(i=1; i<=iter; i++) {
          _cs=c_of[i]
          rm_insertion(_cs)
          gap_length=s_of[i+1] - s_of[i] - gsub(/[MIDSV]/, "", _cs)
          cs=cs c_of[i] padD(gap_length)
        }
        cs=cs c_of[iter+1]
        return cs
      }

    {
      num_of_alignment[$1]++
      start_of[$1]=start_of[$1]","$3
      allele=$4
      reflen=$5
      cstag_of[$1]=cstag_of[$1]","$6
    } END {
      for (read_id in num_of_alignment) {
        sub(/^,/, "" ,start_of[read_id])
        sub(/^,/, "" ,cstag_of[read_id])
        split(start_of[read_id], s_of, ",")
        split(cstag_of[read_id], c_of, ",")
        #* normal
        if (num_of_alignment[read_id]==1) {
          cs=c_of[1]
        }
        #* large deletion
        else if (num_of_alignment[read_id]==2) {
          if (allele !~ /(control)|(wt)/) {
            c_of[1] = padD(gsub(/[MIDS]/, "", c_of[1]))
            c_of[2] = padD(gsub(/[MIDS]/, "", c_of[2]))
          }
          cs=csCat(c_of, s_of, 1)
          }
        #* inversion
        else if (num_of_alignment[read_id]==3) {
          c_of[2] = padD(gsub(/[MIDS]/, "", c_of[2]))
          gsub("D", "V", c_of[2])
          cs=csCat(c_of, s_of, 2)
        }
      print read_id, s_of[1], reflen, cs
    }}' |
    insToNum |
    padding |
    spaceTocomma |
    awk -F, 'NF==$3+3' |
    cut -d, -f 1,4- |
    grep -v "^$" |
    sort -t,
)
