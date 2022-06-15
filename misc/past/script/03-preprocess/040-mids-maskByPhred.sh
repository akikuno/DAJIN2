#!/bin/sh

. $(find .DAJIN_temp/ -name general.sh)
mkdir -p .DAJIN_temp/midsmaskByPhred/ /tmp/.DAJIN_temp/midsmaskByPhred

multi_maskByPhred() {
  cmd=". $(find .DAJIN_temp/ -name maskByPhred.sh); maskByPhred"
  find .DAJIN_temp/sam/"${1:-}"*.sam -print0 |
    xargs -0 -I@ basename @ |
    sed "s/.sam$//" |
    while read -r line; do
      sam=.DAJIN_temp/sam/"$line".sam
      mids=.DAJIN_temp/midsfilter/"$line".csv
      output="$(echo $mids | sed "s|/midsfilter/|/midsmaskByPhred/|")"
      echo "$sam" "$mids" |
        sed "s|^|$cmd |" |
        sed "s|$| >$output \&|"
    done |
    awk -v th="${threads:-1}" '
    {if (NR%(th+1) == 0) print "wait"}
    END {print "wait"}1' |
    sh
}

if find /tmp/.DAJIN_temp/midsmaskByPhred/"$control_name"* 1>/dev/null 2>&1; then
  multi_maskByPhred "$sample_name"
  load_control /tmp/.DAJIN_temp/midsmaskByPhred/
else
  multi_maskByPhred
  save_control .DAJIN_temp/midsmaskByPhred/
fi
