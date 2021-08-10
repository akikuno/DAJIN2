#!/bin/sh

outdir=.DAJIN_temp/consensus/fmtcs/
mkdir -p "$outdir" /tmp/"$outdir"

grep "^@" .DAJIN_temp/sam/"$sample_name"_control.sam >"$outdir"/tmp_header
sort .DAJIN_temp/sam/"$sample_name"_control.sam >"$outdir"/tmp_sam

cat .DAJIN_temp/clustering/"$sample_name".csv |
  cut -d, -f1-2 |
  sort -u |
  while read -r line; do
    suffix="$(echo $line | tr , _)"
    grep "^$line" .DAJIN_temp/clustering/"$sample_name".csv |
      cut -d, -f 3 |
      sort |
      join - "$outdir"/tmp_sam |
      tr " " "\t" |
      cat "$outdir"/tmp_header - >"$outdir"/"$sample_name"_"$suffix".sam
  done

multi_csToCSV() {
  cmd=". $(find .DAJIN_temp/ -name csToCSV.sh); csToCSV"
  find "$outdir"/*.sam |
    while read -r line; do
      output="${line%.*}".csv
      echo "$line" |
        sed "s|^|$cmd |" |
        sed "s|$| >$output \&|"
    done |
    awk -v th="${threads:-1}" '
    {if (NR%(th+1) == 0) print "wait"}
    END {print "wait"}1' |
    sh
}

multi_csToCSV

rm .DAJIN_temp/consensus/fmtcs/*sam
rm .DAJIN_temp/consensus/fmtcs/tmp*

# #? DEBUG =========================================

# cut -d, -f1-2 .DAJIN_temp/clustering/"$sample_name".csv | sort | uniq -c
# wc -l "$outdir"/*
# cat .DAJIN_temp/consensus/fmtcs/barcode31_SV_1.csv | cut -d, -f1 | sort >tmp_id
# grep "^@" .DAJIN_temp/sam/"$sample_name"_control.sam >tmp_header
# sort .DAJIN_temp/sam/"$sample_name"_control.sam |
#   join -v 2 tmp_id - |
#   tr " " "\t" |
#   cat tmp_header - |
#   head >tmp.sam

# # ハードクリップがあるリードの扱いはどうする？
# # 001179a7-d376-42d5-9671-dde7bff87c93

# sort .DAJIN_temp/sam/"$sample_name"_control.sam | awk '$6 ~ /H/' | head
# . library/consensus/csToCSV.sh

# cat tmp.sam |
#   fmtSam |
#   matchToSpace |
#   subToNuc |
#   delToD |
#   awk '{$1=$1","}1' |
#   largeDelAndInvToNuc |
#   insToNuc |
#   padding |
#   spaceTocomma |
#   awk -F, 'NF==$3+3' |
#   cut -d, -f 1,4- |
#   grep -v "^$" |
#   sort -t,

# cat tmp.sam |
#   fmtSam |
#   matchToSpace |
#   subToNuc |
#   delToD |
#   awk '{$1=$1","}1' |
#   largeDelAndInvToNuc |
#   cat >tmp
# cat tmp | grep "G +t Sa St St A"

# echo "id 1 2845 G +t Sa St St A" |
#   insToNuc |
#   padding |
#   spaceTocomma
