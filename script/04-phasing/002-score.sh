#!/bin/sh

cat <<EOF >>log_DAJIN.txt
==========================================================
Classify alleles...
==========================================================
EOF

#----------------------------------------------------------
timestamp "MIDS scoring" >>log_DAJIN.txt
#----------------------------------------------------------

mkdir -p .DAJIN_temp/score /tmp/score

multi_midsToscore() {
  cmd='. .DAJIN_temp/library/midsToscore.sh; midsToscore '
  find .DAJIN_temp/mids/"${1:-}"*.csv |
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
  multi_midsToscore "$sample_name"
  load_control /tmp/score
else
  multi_midsToscore
  save_control .DAJIN_temp/score
fi
