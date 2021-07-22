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
    }'
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
    for(i=2;i<=NF;i++) {
      if($i~/^[acgt]/){
        $i=length($i) $(i+1)
        $(i+1)=""
      }
    }}1'
}

midsconv() (
  cat - |
    cat test/midsconv/input-del_to_wt_allele.sam |
    cat test/midsconv/input-ins_largedel_to_wt_allele.sam |
    fmt_sam |
    matchToM |
    subToS |
    delToD |
    insToI |
    awk '{sub("$", ",", $1)}1' |
    #* Large deletion and Inversion -------------------------
    awk -F, '{
      id_of[$1]++
      flag_of[$1]=flag_of[$1]","$2
      start_of[$1]=start_of[$1]","$3
      allele=$4
      reflen=$5
      cstag_of[$1]=cstag_of[$1]","$6
    } END {
    for(id in id_of) {
      sub(/^,/,"",start_of[id])
      sub(/^,/,"",cstag_of[id])
      split(start_of[id], s_of, ",")
      split(cstag_of[id], c_of, ",")
      #* normal
      if(id_of[id]==1) {
        print id, s_of[1], reflen, c_of[1]
      }
      #* deletion
      else if (id_of[id]==2) {
        if (allele !~ /(control)|(wt)/) {
          str1=""
          gsub(" ", "",c_of[1])
          for (i=1;i<=length(c_of[1]);i++) str1=str1 " D "
          c_of[1] = str1
        }
        cs1=c_of[1]
        gsub(" ", "", cs1)
        gsub("[acgt][acgt]*", "", cs1)
        del_len=s_of[2]-s_of[1]-length(cs1)
        str=""
        for(i=1;i<=del_len;i++) str=str " D "
        cs=c_of[1]" "str" "c_of[2]
        print id, s_of[1], reflen, cs
        }
      #* inversion
      else if (id_of[id]==3) {
        if (allele !~ /(control)|(wt)/) {
          str=""
          gsub(" ", "",c_of[1])
          for (i=1;i<=length(c_of[1]);i++) str=str " D "
          c_of[1] = str
        }
        str1=""; str2=""
        cs1=c_of[1]; gsub(" ", "", cs1); gsub("[acgt][acgt]*", "", cs1)
        cs2=c_of[2]; gsub(" ", "", cs2); gsub("[acgt][acgt]*", "", cs2)
        del_len1=s_of[2]-s_of[1]-length(cs1)
        del_len2=s_of[3]-s_of[2]-length(cs2)
        for(i=1;i<=del_len1;i++) str1=str1 " D "
        for(i=1;i<=del_len2;i++) str2=str2 " D "
        cs=c_of[1]" "str1" "c_of[2]" "str2" "c_of[3]
        print id, s_of[1], reflen, cs
      }
    }}' |
    sed "s/  */ /g" |
    sed "s/  */ /g" |
    #* Padding ----------------------------------------------
    awk '{
      start=$2-1; end=$3; len=NF-3+start
      pad_start=""; pad_end=""
      for(i=1;i<=start;i++) pad_start=pad_start "= "
      $4=pad_start $4
      for(i=len;i<end;i++) pad_end=pad_end " ="
      $NF=$NF pad_end
    }1' |
    sed "s/  */ /g" |
    awk 'NF-3 == $3' |
    #* Convert to CSV ---------------------------------------
    sed "s/ /,/g" |
    sed "s/,$//" |
    cut -d, -f 1,4- |
    grep -v "^$" |
    sort -t,
)
