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

load_control() {
  find "$1" -type f |
    while read -r line; do
      output=${line#/tmp/}
      output=${output%.gz}
      gzip -dc "$line" >.DAJIN_temp/"$output"
    done
}

save_control() {
  find "$1" -type f |
    grep "$control_name" |
    while read -r line; do
      output=${line#\.DAJIN_temp/}.gz
      gzip -c "$line" >/tmp/"$output"
    done
}

if find /tmp/scalar/"$control_name"* 1>/dev/null 2>&1; then
  multi_midsToscore "$sample_name"
  load_control /tmp/scalar
else
  multi_midsToscore
  save_control .DAJIN_temp/scalar
fi
