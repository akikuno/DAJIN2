#!/bin/sh

#----------------------------------------------------------
timestamp "MIDS encoding" >>log_DAJIN.txt
#----------------------------------------------------------

mkdir -p .DAJIN_temp/midsmask /tmp/DAJIN/midsmask

multi_maskMIDS() {
  cmd='. .DAJIN_temp/library/maskMIDS.sh; maskMIDS '
  find .DAJIN_temp/sam/"${1:-}"*.sam -print0 |
    xargs -0 -I@ basename @ |
    sed "s/.sam$//" |
    while read -r line; do
      sam=.DAJIN_temp/sam/"$line".sam
      mids=.DAJIN_temp/mids/"$line".csv
      output="$(echo $mids | sed "s|/mids/|/midsmask/|")"
      echo "$sam" "$mids" |
        sed "s|^|$cmd|" |
        sed "s|$| >$output \&|"
    done |
    awk -v th="${threads:-1}" '
    {if (NR%(th+1) == 0) print "wait"}
    END {print "wait"}1' |
    sh
}

if find /tmp/DAJIN/midsmask/"$control_name"* 1>/dev/null 2>&1; then
  multi_maskMIDS "$sample_name"
  load_control /tmp/DAJIN/midsmask/
else
  multi_maskMIDS
  save_control .DAJIN_temp/midsmask/
fi
