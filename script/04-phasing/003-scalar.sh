#!/bin/sh

mkdir -p .DAJIN_temp/scalar /tmp/scalar

multiscoreToScalar() {
  cmd='. .DAJIN_temp/library/scoreToScalar.sh; scoreToScalar '
  find .DAJIN_temp/score/"$1"* |
    grep -e "$sample_name" -e "$control_name"_wt -e "$control_name"_control |
    while read -r line; do
      output="$(echo $line | sed "s|/score/|/scalar/|")"
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
  multiscoreToScalar "$sample_name"
  load_control /tmp/scalar
else
  multiscoreToScalar
  save_control .DAJIN_temp/scalar
fi
