#!/bin/sh

#----------------------------------------------------------
timestamp "MIDS encoding" >>log_DAJIN.txt
#----------------------------------------------------------

mkdir -p .DAJIN_temp/mids /tmp/mids

multi_samTomids() {
  cmd='. .DAJIN_temp/library/samTomids.sh; samTomids '
  find .DAJIN_temp/sam/"$1"*.sam |
    while read -r line; do
      output="${line%.*}".csv
      output="$(echo $output | sed "s|/sam/|/mids/|")"
      echo "$line" |
        sed "s|^|$cmd|" |
        sed "s|$| >$output \&|"
    done |
    awk -v th="${threads:-1}" '
    {if (NR%(th+1) == 0) print "wait"}
    END {print "wait"}1' |
    sh
}

if find /tmp/mids/"$control_name"* 1>/dev/null 2>&1; then
  find /tmp/mids/* |
    while read -r line; do
      gzip -dc "$line" >.DAJIN_temp/mids/"$(basename ${line%.gz})"
    done
else
  multi_samTomids "$control_name"
  find .DAJIN_temp/mids/"$control_name"* |
    while read -r line; do
      gzip -c "$line" >/tmp/mids/"$(basename $line)".gz
    done
fi

multi_samTomids "$sample_name"
