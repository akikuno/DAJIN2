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
