#!/bin/sh

error_exit() {
  echo "ERROR: $1" >&2
  echo "ERROR: $1" >>log_DAJIN.txt
  [ "${TEST:-}" ] || rm -rf .DAJIN_temp/ 2>/dev/null
  exit 1
}

terminate() {
  trap '' TERM
  kill -TERM 0
  [ "${TEST:-}" ] || rm -rf .DAJIN_temp/ 2>/dev/null
  exit "$1"
}
trap "terminate 130" INT
trap "terminate 143" TERM

timestamp() {
  echo "$(date +'%Y-%m-%d %H:%M:%S') | $*"
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

open() {
  if file "$1" | grep -q gzip; then
    gzip -dc "$1"
  else
    cat "$1"
  fi
}
