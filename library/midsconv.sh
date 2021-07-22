#!/bin/sh
# shellcheck disable=SC1091,SC2120

################################################################################
# input: SAM file using `minimap2 -ax splice --cs=long`
# output: Match, Insertion, Deletion, Substitution, and "= (null)"
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

insToI() {
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

midsconv() (
  cat - |
    # cat test/midsconv/input-del_to_wt_allele.sam |
    # cat test/midsconv/test_inv.sam |
    fmt_sam |
    matchToM |
    subToS |
    delToD |
    awk '{$1=$1","}1' |
    #* Large deletion and Inversion -------------------------
    awk -F, '
      function pad_D (_length) {
        str=""
        for (i=1;i<=_length;i++) str=str " D "
        return str
      }
      function rm_insertion(_string) {
        gsub("[acgt][acgt]*", "", _string)
        return _string
      } {
      num_of_alignment[$1]++
      start_of[$1]=start_of[$1]","$3
      allele=$4
      reflen=$5
      cstag_of[$1]=cstag_of[$1]","$6
    } END {
    for (id in num_of_alignment) {
      sub(/^,/, "" ,start_of[id])
      sub(/^,/, "" ,cstag_of[id])
      split(start_of[id], s_of, ",")
      split(cstag_of[id], c_of, ",")
      #* normal
      if (num_of_alignment[id]==1) {
        print id, s_of[1], reflen, c_of[1]
      }
      #* large deletion
      else if (num_of_alignment[id]==2) {
        #* controlアレル以外のアレルの構造多型はすべて異常としたいため.
        if (allele !~ /(control)|(wt)/) {
          c_of[1] = pad_D(length(c_of[1]))
        }
        cs1=c_of[1]; rm_insertion(cs1)
        D=pad_D(s_of[2] - s_of[1] - length(cs1))
        cs=c_of[1] " " D " " c_of[2]
        print id, s_of[1], reflen, cs
        }
      #* inversion
      else if (num_of_alignment[id]==3) {
        if (allele !~ /(control)|(wt)/) {
          c_of[1] = pad_D(c_of[1])
        }
        # for(i=i; i<=2; i++) {
        #   cs=c_of[i]
        #   rm_insertion(cs)
        #   array[i]=pad_D(s_of[i+1]-s_of[i]-length(cs))
        # }
        cs1=c_of[1]; rm_insertion(cs1)
        cs2=c_of[2]; rm_insertion(cs2)
        D1 = pad_D(s_of[2]-s_of[1]-length(cs1))
        D2 = pad_D(s_of[3]-s_of[2]-length(cs2))
        cs=c_of[1] " " D1 " " c_of[2] " " D2 " " c_of[3]
        print id, s_of[1], reflen, cs
      }
    }}' |
    insToI |
    padding |
    spaceTocomma |
    cut -d, -f 1,4- |
    grep -v "^$" |
    sort -t,
)
