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

matchToSpace() {
  awk '{
    for(i=2;i<=NF;i++)
      gsub(/[ACGT]/, "& ", $i)
  }1'
}

subToNuc() {
  awk '{
    for (i=2;i<=NF;i++) {
      if ($i ~ /^\*/) {
        gsub(/\*[acgt]/, "S", $i)
      }
  }}1'
}

delToD() {
  awk '{
    for(i=2;i<=NF;i++) {
      if($i ~ /^-/) {
        str=""
        for(j=1; j<=int(length($i)-1); j++) str=str " D "
        $i=str
      }
    }}1'
}

largeDelAndInvToNuc() {
  awk -F, '
      function padD(iter,    i,str) {
        for (i=1; i<=iter; i++) str=str " D "
        return str
      }

      function ins_rm(string) {
        gsub("[acgt][acgt]*", "", string)
        return string
      }

      function csCat(c_of, s_of, iter,    i,cs) {
        for(i=1; i<=iter; i++) {
          _cs=c_of[i]
          ins_rm(_cs)
          gap_length=s_of[i+1] - s_of[i] - gsub(/[ACGTacgt]/, "", _cs)
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
          cs=csCat(c_of, s_of, 1)
          }
        #* inversion
        else if (num_of_alignment[read_id]==3) {
          c_of[2] = tolower(c_of[2])
          cs=csCat(c_of, s_of, 2)
        }
      print read_id, s_of[1], reflen, cs
    }}'
}

insToNuc() {
  awk '{
    for (i=4; i<=NF; i++) {
      if ($i ~ /^\+/){
        sub(/\+/, "", $i)
        $i = "I" $i substr($(i+1), 1, 1)
        sub(/[ACGT]/, "", $(i+1))
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
    for (i=1;i<=start;i++)
      pad_start = pad_start " = "
    for (i=len;i<end;i++)
      pad_end = pad_end " = "
    $4 = pad_start $4
    $NF = $NF pad_end
  }1' |
    sed "s/  */ /g"
}

spaceTocomma() {
  sed -e "s/  */,/g" -e "s/,$//"
}

csToCSV() (
  if [ -p /dev/stdin ] && [ "$#" -eq 0 ]; then
    cat -
  elif [ -r "$1" ]; then
    cat "$1"
  else
    echo "$*"
  fi |
    # cat test/csToCSV/que.sam |
    fmtSam |
    matchToSpace |
    subToSpace |
    delToSpace |
    awk '{$1=$1","}1' |
    #* Large deletion and Inversion -------------------------
    awk -F, '
      function padD(iter,    i,str) {
        for (i=1; i<=iter; i++) str=str "D"
        return str
      }

      function ins_rm(string) {
        gsub("[acgt][acgt]*", "", string)
        return string
      }

      function csCat(c_of, s_of, iter,    i,cs) {
        for(i=1; i<=iter; i++) {
          _cs=c_of[i]
          ins_rm(_cs)
          gap_length=s_of[i+1] - s_of[i] - gsub(/[ACGTacgt]/, "", _cs)
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
          cs=csCat(c_of, s_of, 1)
          }
        #* inversion
        else if (num_of_alignment[read_id]==3) {
          c_of[2] = tolower(c_of[2])
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
