#!/bin/sh

mkdir -p .DAJIN_temp/scalar /tmp/.DAJIN_temp/scalar

multiscoreToScalar() {
  cmd=". $(find .DAJIN_temp/ -name scoreToScalar.sh); scoreToScalar"
  find .DAJIN_temp/score/"${1:-}"* |
    grep -e "$sample_name" -e "$control_name"_wt -e "$control_name"_control |
    while read -r line; do
      output="$(echo $line | sed "s|/score/|/scalar/|")"
      echo "$line" |
        sed "s|^|$cmd |" |
        sed "s|$| >$output \&|"
    done |
    awk -v th="${threads:-1}" '
    {if (NR%(th+1) == 0) print "wait"}
    END {print "wait"}1' |
    sh
}

if find /tmp/.DAJIN_temp/scalar/"$control_name"* 1>/dev/null 2>&1; then
  multiscoreToScalar "$sample_name"
  load_control /tmp/.DAJIN_temp/scalar
else
  multiscoreToScalar
  save_control .DAJIN_temp/scalar
fi
