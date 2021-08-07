#!/bin/sh
# shellcheck disable=SC1090,SC1091,SC2002,SC2030,SC2086,SC2097,SC2098,SC2046

# DAJIN consensus -t [threads] -a [fasta] -c [fasta|fastq|bam|sam] -s [fasta|fastq|bam|sam|CRAM] -f [html|fasta] -s [strand]
# -c|--control
# -s|--sample
# input: SAM wich CS-tag
# output: FASTA/HTML consensus sequence

################################################################################
# Initialization
################################################################################

set -eu

error_exit() {
  echo "ERROR: $1" >&2
  echo "ERROR: $1" >>log_DAJIN.txt
  # rm -rf .DAJIN_temp/ 2>/dev/null
  exit 1
}

terminate() {
  trap '' TERM
  kill -TERM 0
  # rm -rf .DAJIN_temp/ 2>/dev/null
  exit "$1"
}
trap "terminate 130" INT
trap "terminate 143" TERM

mkdir -p .DAJIN_temp/consensus

#----------------------------------------------------------
# Preprocessing
# output: id, location, ref nuc, sample nuc, sample mids
#----------------------------------------------------------

. .DAJIN_temp/library/sepconv.sh
### TEMPORAL
: "${format:=html}"
ref_fasta=$1
cont_consensus=$2
sample_consensus=$3
: "${strand:=-}"=$4
threads=$(getconf _NPROCESSORS_ONLN | awk '{print $0-2}')
sample_consensus_name="$(echo ${sample_consensus##*/} | cut -d_ -f 1)"
###

################################################################################
# コントロールとサンプルでMIDS頻度を引き算をします.
################################################################################

: （もしなかったら）コントロールのMIDS頻度を作ります && {
  if ! [ -f .DAJIN_temp/consensus/midsfreq_control.csv ]; then
    if samtools view -H "$cont_consensus" >/dev/null 2>&1; then
      if samtools view "$cont_consensus" | head -n 100 | grep -q "cs:Z:="; then
        cp "$cont_consensus" .DAJIN_temp/consensus/tmp_control.sam
      else
        samtools view "$cont_consensus" |
          awk '$2==0 || $2==16 {print ">"$1"\n"$10}' |
          minimap2 -t "${threads}" -ax map-ont "${ref_fasta}" - --cs=long 2>/dev/null |
          cat >.DAJIN_temp/consensus/tmp_control.sam
      fi
    else
      minimap2 -t "${threads}" -ax map-ont "${ref_fasta}" "$cont_consensus" --cs=long 2>/dev/null |
        cat >.DAJIN_temp/consensus/tmp_control.sam
    fi
    [ -s .DAJIN_temp/consensus/tmp_control.sam ] || error_exit "control file was unrecognized format"

    . .DAJIN_temp/library/midsconv.sh
    midsconv .DAJIN_temp/consensus/tmp_control.sam >.DAJIN_temp/consensus/tmp_mids_control.csv
    [ -s .DAJIN_temp/consensus/tmp_mids_control.csv ] || error_exit "control file was unrecognized format"

    . .DAJIN_temp/library/calc_mids_freq.sh
    cat .DAJIN_temp/consensus/tmp_mids_control.csv |
      cut -d, -f 2- |
      sed "s/^/control,/" |
      sed "s/,=/,M/g" |
      sed "s/,[0-9][0-9]*[MDS]/,I/g" |
      calc_mids_freq |
      sed "s/,/@/" |
      sort -u -t, |
      cat >.DAJIN_temp/consensus/midsfreq_control.csv
    rm .DAJIN_temp/consensus/tmp*
  fi
}

#----------------------------------------------------------
# サンプルのMIDS頻度を用意します
#----------------------------------------------------------

mids_converted=$(
  find .DAJIN_temp/mids/"$sample_consensus_name"*.csv |
    grep -e _wt.csv -e _control.csv
)

if [ "${mids_converted:-}" ]; then
  : MIDS変換済みのファイルがあればそれを使います. && {
    grep -v "^@" "$sample_consensus" |
      cut -f 1 |
      sort -u |
      join -t, "$mids_converted" - >.DAJIN_temp/consensus/tmp_mids_sample.csv
  }
else
  : MIDS変換済みのファイルがなければ生成します && {
    if samtools view -H "$sample_consensus" >/dev/null 2>&1; then
      if samtools view "$sample_consensus" | head -n 100 | grep -q "cs:Z:="; then
        cp "$sample_consensus" .DAJIN_temp/consensus/tmp_sample.sam
      else
        samtools view "$sample_consensus" |
          awk '$2==0 || $2==16 {print ">"$1"\n"$10}' |
          minimap2 -t "${threads}" -ax map-ont "${ref_fasta}" - --cs=long 2>/dev/null |
          cat >.DAJIN_temp/consensus/tmp_sample.sam
      fi
    else
      minimap2 -t "${threads}" -ax map-ont "${ref_fasta}" "$sample_consensus" --cs=long 2>/dev/null |
        cat >.DAJIN_temp/consensus/tmp_sample.sam
    fi
    [ -s .DAJIN_temp/consensus/tmp_sample.sam ] || error_exit "sample file was unrecognized format"
    . .DAJIN_temp/library/midsconv.sh
    midsconv .DAJIN_temp/consensus/tmp_sample.sam >.DAJIN_temp/consensus/tmp_mids_sample.csv
  }
fi

. .DAJIN_temp/library/calc_mids_freq.sh
cat .DAJIN_temp/consensus/tmp_mids_sample.csv |
  cut -d, -f 2- |
  sed "s/^/control,/" |
  sed "s/,=/,M/g" |
  sed "s/,[0-9][0-9]*[MDS]/,I/g" |
  calc_mids_freq |
  sed "s/,/@/" |
  sort -u -t, |
  cat >.DAJIN_temp/consensus/midsfreq_sample.csv

#----------------------------------------------------------
# サンプルからコントロールのMIDS頻度を引き算して, シークエンスエラーを考慮したコンセンサス配列を用意します
#----------------------------------------------------------

cat .DAJIN_temp/consensus/midsfreq_sample.csv |
  join -t, - .DAJIN_temp/consensus/midsfreq_control.csv |
  awk -F, 'function RELU(v) {return v < 0 ? 0 : v}
    BEGIN {OFS=","} {
    subI=RELU($3-$7)
    subD=RELU($4-$8)
    subS=RELU($5-$9)
    subM=RELU(1-subI-subD-subS)
    # print subM,subI,subD,subS
    if ($3>0.9) MIDS="I"
    else if ($4>0.9) MIDS="D"
    else if ($5>0.9) MIDS="S"
    #* sequence error
    else if (subI>subM && subI>subD && subI>subS) MIDS="I"
    else if (subD>subM && subD>subI && subD>subS) MIDS="D"
    else if (subS>subM && subS>subI && subS>subD) MIDS="S"
    else MIDS="M"
    print $0,MIDS
  }' |
  tr "@" "," |
  awk -F, '{print $2,$NF}' |
  sort -n |
  awk '{printf $NF","}' |
  sed "s/,$//" |
  grep ^ >.DAJIN_temp/consensus/tmp_mids_consensus.csv
[ -s .DAJIN_temp/consensus/tmp_mids_consensus.csv ] || error_exit "tmp_mids_consensus.csv was empty"

#----------------------------------------------------------
# Obtain base information corresponding to MIDS
# サンプルの塩基情報を入手します.
#----------------------------------------------------------

rm .DAJIN_temp/consensus/tmp_split_* 2>/dev/null || :
split_num=$(awk -v th="$((threads - 1))" 'END{print int(NR/th)}' "$sample_consensus")
split -l "$split_num" "$sample_consensus" .DAJIN_temp/consensus/tmp_split_

grep "^@" "$sample_consensus" >.DAJIN_temp/consensus/tmp_header
find .DAJIN_temp/consensus/tmp_split_* |
  while read -r line; do
    cat .DAJIN_temp/consensus/tmp_header "$line" >"$line"_sam
  done

cmd='. .DAJIN_temp/library/sepconv.sh && sepconv'
find .DAJIN_temp/consensus/tmp_split_*_sam |
  awk -v cmd="$cmd" '{
    sample1=$0; output=$0
    sub("sam$","csv",output)
    print cmd, sample1, ">", output"_sepconv", "&"
    } END {print "wait"}' |
  sh

cat .DAJIN_temp/consensus/tmp_split*_sepconv >.DAJIN_temp/consensus/tmp_sepconv.csv

#----------------------------------------------------------
# サンプルの塩基情報とMIDSコンセンサスをもとに, IDS部分の塩基情報を得ます.
#----------------------------------------------------------

. .DAJIN_temp/library/consensus.sh

cut -d, -f2- .DAJIN_temp/consensus/tmp_sepconv.csv |
  awk -F, '
    BEGIN {OFS=","
      getline mids < ".DAJIN_temp/consensus/tmp_mids_consensus.csv"
      n=split(mids, mids_of, ",") } {
    for(i=1;i<=n;i++) {
    #* MIDSがInsの場合, sepconvもIns以外を消去します.
      if (mids_of[i] == "I") {
        if ($i !~ /^[acgt][acgt]*[ACGT]/) $i=""
      }
    }}1' |
  # awk -F, '
  #   BEGIN {OFS=","
  #     getline mids < ".DAJIN_temp/consensus/tmp_mids_consensus.csv"
  #     n=split(mids, mids_of, ",") } {
  #   for(i=1;i<=n;i++) {
  #   #* Insertionの塩基数を求めます.
  #     if (mids_of[i] ~ /^[0-9]+[A-Z]/) {
  #       sub(/[A-Z]/, "", mids_of[i])
  #       mids_of[i]++
  #     }
  #   #* Insertion以外は塩基数は1になります.
  #     else {
  #       mids_of[i]=1
  #     }
  #   #* Insertionの塩基数がtmp_mids_consensus.csvとtmp_sepconv.csvで違う場合, シークエンスエラーとみなして無視します
  #     if (length($i) != mids_of[i]) $i=""
  #   }}1' |
  consensus |
  awk '{print NR,$0}' |
  sort -t " " >.DAJIN_temp/consensus/tmp_mids_sep.ssv

#----------------------------------------------------------
# 参照配列と結合します.
#----------------------------------------------------------

sed 1d "$ref_fasta" |
  awk '{n=split($0,array,""); for(i=1;i<=n;i++) print i,array[i]}' |
  sort -t" " >.DAJIN_temp/consensus/tmp_ref.ssv

cat .DAJIN_temp/consensus/tmp_mids_consensus.csv |
  tr "," "\n" |
  awk '{print NR,$0}' |
  sort -t " " |
  join .DAJIN_temp/consensus/tmp_mids_sep.ssv - |
  join .DAJIN_temp/consensus/tmp_ref.ssv - |
  #* insertion
  sed "s/[0-9][0-9]*M$/I/" |
  #* the same frequency
  sed "s/= S/= M/" |
  #* padding
  tr "=" "D" |
  tr "acgt" "ACGT" |
  sort -n |
  #* inversion
  awk '
    $3=="V" {
      row[NR]=$0
      if(min_nr=="") min_nr=NR
      if(max_nr<NR)  max_nr=NR
      next
    } {print}
    END {
      nr=min_nr
      for(i=max_nr; i>=min_nr;i--){
        $0=row[i]
        $1=nr
        if($2=="A") $2="T"
        else if($2=="C") $2="G"
        else if($2=="G") $2="C"
        else if($2=="T") $2="A"
        nr++
        print
      }
    }' |
  grep -v "^$" |
  tr " " "," |
  sort -t, -n >.DAJIN_temp/consensus/tmp_loc_ref_sample_mids.csv

: アンチセンスの場合は逆相補鎖とします && {
  if [ _"$strand" = "_-" ]; then
    sort -t, -nr .DAJIN_temp/consensus/tmp_loc_ref_sample_mids.csv |
      tr "ACGT" "TGCA" |
      awk -F, 'BEGIN {OFS=","}
        $NF=="I" {
          rev=""
          n=split($3,array,"")
          for(i=n;i>0;i--) rev=rev array[i]
          $3=rev
          print
          next}1' |
      awk -F, 'BEGIN{OFS=","}{$1=NR}1' >.DAJIN_temp/consensus/tmp
    mv .DAJIN_temp/consensus/tmp .DAJIN_temp/consensus/tmp_loc_ref_sample_mids.csv
  fi
}

fasta_id="$(basename "${sample_consensus%.sam}")"

: FASTA REPORTS && {
  cat .DAJIN_temp/consensus/tmp_loc_ref_sample_mids.csv |
    grep -v "D" |
    awk -F, '{
      if ($4 ~/[SI]/) $2=$3
      print $2
    }' |
    tr -d "\n" |
    awk -v id="$fasta_id" 'BEGIN {print ">"id}1' |
    grep ^ >"${sample_consensus%.sam}".fasta
}

: HTML REPORTS && {
  cat .DAJIN_temp/consensus/tmp_loc_ref_sample_mids.csv |
    awk -F, 'NF==4' |
    awk -F, '{
      if ($4=="S") {$2="<span class=\"Sub\">"$3"</span>"}
      else if ($4=="I") {
        INS=substr($3, 1, length($3)-1)
        Match=substr($3, length($3), 1)
        $2="<span class=\"Ins\">" INS "</span>" Match
        }
      else if ($4=="D") {$2="<span class=\"Del\">"$2"</span>"}
      else {$2=$2}
      if ($3=="V") {$2="<span class=\"Inv\">"$2"</span>"}
      print $2
    }' |
    tr -d "\n" |
    sed "s/^/>${fasta_id}@@@/" |
    cat .DAJIN_temp/utils/consensus_head.html - .DAJIN_temp/utils/consensus_tail.html |
    awk '{gsub("@@@","\n")}1' |
    grep ^ >"${sample_consensus%.sam}".html
}

mv .DAJIN_temp/consensus/tmp_loc_ref_sample_mids.csv "${sample_consensus%.sam}"_loc_ref_sample_mids.csv

# #TODO LOXPでの検討が必要. そもそも必要かな？
# : Rename mutant allele && {
#   find .DAJIN_temp/consensus/"$sample_consensus_name"_*.fa |
#     grep -v -e "SV.fa" |
#     while read -r line; do
#       alleletype=$(echo "${line#*${sample_consensus_name}_}" | sed "s/.fa//")
#       sed 1d "$line" >.DAJIN_temp/consensus/tmp_sample
#       sed 1d .DAJIN_temp/input/"$alleletype" >.DAJIN_temp/consensus/tmp_control
#       if ! diff .DAJIN_temp/consensus/tmp_sample .DAJIN_temp/consensus/tmp_control 1>/dev/null 2>&1; then
#         bf_header="$(head -n 1 "$line")"
#         af_header="$(echo ">${line##*/}" | sed "s/.fa$//")"
#         : _mutated_を追加する: && {
#           rename="${line%_*}"
#           rename="${rename%_mutant*}"_mutant_"$alleletype"
#         }
#         if [ _"$line" != _"$rename" ]; then
#           # FASTA -----------------
#           sed "s/${bf_header}/${af_header}/" "$line" |
#             grep ^ >"$rename"
#           # HTML -----------------
#           sed "s/${bf_header}/${af_header}/" "${line%.fa}.html" |
#             grep ^ >"${rename%.fa}".html
#           # SAM -----------------
#           mv "${line%.fa}".sam "${rename%.fa}".sam 2>/dev/null || :
#           # Clearn --------------
#           rm "${line}" "${line%.fa}".html
#         fi
#       fi
#     done
# }
