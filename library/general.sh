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
