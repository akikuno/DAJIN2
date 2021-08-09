#!/bin/sh
# shellcheck disable=SC1091,SC2120

################################################################################
# input: SAM file using `minimap2 -ax splice --cs=long`
# output: DNA sequence with MIDS separation
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

grep -e "^@" -e "cs:Z" |
  #* Get id,start,length, converted CStag (eg: A+cg*cG-c => A c@cG@c)
  awk '
      /^@SQ/ && /LN:/ {
        for(i=1;i<=NF;i++) if($i ~ /^LN:/) len=$i
        sub(/.*LN:/,"", len); next}
      /cs:Z:/ {
        id=$1;start=$4; cstag=$(NF-1)
        sub("cs:Z:","",cstag)
        gsub("=","",cstag)
        gsub(/\-/,"@@@",cstag)
        gsub(/\*[acgt]/,"@@@",cstag)
        gsub(/\+/," ",cstag)
        gsub("~"," ~",cstag)
      print id","start","len" "cstag}' |
  awk '{sub("$",",",$1)}1' |
  sed "s/  */ /g" |
  sed "s/, /,/g" |
  sort -t, -k 2,2n |
  #* Large deletion and Inversion -------------------------
  awk -F, '{
      id_of[$1]++
      start_of[$1]=start_of[$1]","$2
      reflen=$3
      cstag_of[$1]=cstag_of[$1]","$4
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
        str=""
        cs1=c_of[1]; gsub(" ", "", cs1); gsub("@@@", "", cs1)
        del_len=s_of[2]-s_of[1]-length(cs1)+2
        for(i=1;i<=del_len;i++) str=str "D"
        cs=c_of[1]" "str" "c_of[2]
        print id, s_of[1], reflen, cs
        }
      #* inversion
      else if (id_of[id]==3) {
        str1=""; str2=""; cs2_inv=""
        cs1=c_of[1]; gsub(" ", "", cs1); gsub("@@@", "", cs1)
        cs2=c_of[2]; gsub(" ", "", cs2); gsub("@@@", "", cs2)
        #* annotate inv as V
        gsub(/[ACGT]/, "V", c_of[2])
        del_len1=s_of[2]-s_of[1]-length(cs1)+2
        del_len2=s_of[3]-s_of[2]-length(cs2)+2
        for(i=1;i<=del_len1;i++) str1=str1 "D"
        for(i=1;i<=del_len2;i++) str2=str2 "D"
        cs=c_of[1]" "str1" "c_of[2]" "str2" "c_of[3]
        print id, s_of[1], reflen, cs
      }
    }}' |
  awk '{$1=$1","; $2=$2","}1' |
  sed "s/  */ /g" |
  sed "s/, /,/g" |
  #* CSV format -------------------------------------------
  awk 'BEGIN{OFS=","}{
      for(i=2;i<=NF;i++){
      seq=""
      match($i,/^[acgt]+/)
      #* Insertion (eg: acgGAA => acgG,A,A; ac@gA => acg,A)
      if(RSTART==1) {
        ins_len=RLENGTH+1
        rest_len=RLENGTH+2
        if($i~/^[acgt][acgt]*@@@/) {
          ins_len=RLENGTH+4
          rest_len=RLENGTH+5
        }
        ins=substr($i,1,ins_len)","
        rest=substr($i,rest_len)
        split(rest,array,"")
        for(j in array) seq=seq array[j]","
        $i=ins seq
      } else {
      #* Other (eg: AcaG => A,c,a,G)
        split($i,array,"")
        for(j in array) seq=seq array[j]","
        $i=seq
      }
    }}1' |
  sed "s/@,@,@//g" |
  sed "s/@@@//g" |
  sed "s/,,/,/g" |
  sed "s/,$//" |
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
