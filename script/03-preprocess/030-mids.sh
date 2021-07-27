#!/bin/sh

#----------------------------------------------------------
echo "$(date +'%Y-%m-%d %H:%M:%S') MIDSV encoding" >>log_DAJIN.txt
#----------------------------------------------------------

mkdir -p .DAJIN_temp/midsv
cmd='. .DAJIN_temp/library/samTomids.sh; samTomids '

if find /tmp/"$control_name"*.csv 1>/dev/null 2>&1; then
  cp /tmp/"$control_name"*.csv .DAJIN_temp/midsv
else
  find .DAJIN_temp/sam/"$control_name"*.sam |
    while read -r line; do
      output="${line%.*}".csv
      echo "$line" |
        sed "s|^|$cmd|" |
        sed "s|$| >$output \&|"
    done |
    awk -v th="$threads" '
    {if (NR%(th+1) == 0) print "wait"}
    END {print "wait"}1' |
    sh
  cp .DAJIN_temp/sam/"$control_name"*.csv /tmp/
fi

find .DAJIN_temp/sam/"$sample_name"*.sam |
  while read -r line; do
    output="${line%.*}".csv
    echo "$line" |
      sed "s|^|$cmd|" |
      sed "s|$| >$output \&|"
  done |
  awk -v th="$threads" '
    {if (NR%(th+1) == 0) print "wait"}
    END {print "wait"}1' |
  sh

mv .DAJIN_temp/sam/*.csv .DAJIN_temp/midsv
