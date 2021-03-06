#!/bin/sh

cat <<EOF >>log_DAJIN.txt
==========================================================
Classify alleles...
==========================================================
EOF

#----------------------------------------------------------
timestamp "MIDS scoring" log_DAJIN.txt
#----------------------------------------------------------

mkdir -p .DAJIN_temp/score /tmp/.DAJIN_temp/score

multi_midsToScore() {
  cmd=". $(find .DAJIN_temp/ -name midsToScore.sh); midsToScore"
  find .DAJIN_temp/midsmaskBySubtract/"${1:-}"*.csv |
    while read -r line; do
      output="${line%.*}".csv
      output="$(echo $output | sed "s|/midsmaskBySubtract/|/score/|")"
      echo "$line" |
        sed "s|^|$cmd |" |
        sed "s|$| >$output \&|"
    done |
    awk -v th="${threads:-1}" '
    {if (NR%(th+1) == 0) print "wait"}
    END {print "wait"}1' |
    sh
}

if find /tmp/.DAJIN_temp/score/"$control_name"* 1>/dev/null 2>&1; then
  multi_midsToScore "$sample_name"
  load_control /tmp/.DAJIN_temp/score
else
  multi_midsToScore
  save_control .DAJIN_temp/score
fi
