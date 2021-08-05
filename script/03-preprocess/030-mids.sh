#!/bin/sh

#----------------------------------------------------------
timestamp "MIDS encoding" >>log_DAJIN.txt
#----------------------------------------------------------

mkdir -p .DAJIN_temp/mids /tmp/DAJIN/mids

multi_samToMIDS() {
  cmd='. .DAJIN_temp/library/samToMIDS.sh; samToMIDS '
  find .DAJIN_temp/sam/"${1:-}"*.sam |
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

if find /tmp/DAJIN/mids/"$control_name"* 1>/dev/null 2>&1; then
  multi_samToMIDS "$sample_name"
  load_control /tmp/DAJIN/mids
else
  multi_samToMIDS
  save_control .DAJIN_temp/mids
fi
