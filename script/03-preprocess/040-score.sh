#!/bin/sh

#----------------------------------------------------------
timestamp "MIDS scoring" >>log_DAJIN.txt
#----------------------------------------------------------

mkdir -p .DAJIN_temp/score /tmp/score

multi_midsToscore() {
  cmd='. .DAJIN_temp/library/midsToscore.sh; midsToscore '
  find .DAJIN_temp/mids/"$1"*.csv |
    while read -r line; do
      output="${line%.*}".csv
      output="$(echo $output | sed "s|/mids/|/score/|")"
      echo "$line" |
        sed "s|^|$cmd|" |
        sed "s|$| >$output \&|"
    done |
    awk -v th="${threads:-1}" '
    {if (NR%(th+1) == 0) print "wait"}
    END {print "wait"}1' |
    sh
}

if find /tmp/score/"$control_name"* 1>/dev/null 2>&1; then
  find /tmp/score/* |
    while read -r line; do
      gzip -dc "$line" >.DAJIN_temp/score/"$(basename ${line%.gz})"
    done
else
  multi_midsToscore "$control_name"
  find .DAJIN_temp/score/"$control_name"* |
    while read -r line; do
      gzip -c "$line" >/tmp/score/"$(basename $line)".gz
    done
fi

multi_midsToscore "$sample_name"
