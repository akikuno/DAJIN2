#!/bin/sh

set -ue
umask 0022
export LC_ALL=C

[ "${BASH_VERSION:-}" ] && set -o posix -o pipefail
[ "${ZSH_VERSION:-}" ] && setopt shwordsplit interactivecomments

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
  echo "$(date +'%Y-%m-%d %H:%M:%S') $*"
}

if type wget >/dev/null 2>&1; then
  CMD_CHECK='wget -q -O - --spider --tries=2 --wait=1 --timeout=5'
  CMD_GET='wget -q -O -'
elif type curl >/dev/null 2>&1; then
  CMD_CHECK='curl --retry 2 --retry-delay 1 -s -o /dev/null -w "%{http_code}"'
  CMD_GET='curl -s'
else
  error_exit 'No HTTP-GET/POST command found.'
fi

timestamp "DAJIN $*" >log_DAJIN.txt

find ./ -type f |
  grep -e lib -e src -e doc -e utils |
  grep -v -e "init" -e "past" -e ".DAJIN_temp" |
  cut -d / -f 1-2 |
  sort -u |
  sed "s|^./|.DAJIN_temp/|" |
  xargs mkdir -p
