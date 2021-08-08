#!/bin/sh

set -eu
umask 0022
export LC_ALL=C

[ "${BASH_VERSION:-}" ] && set -o posix -o pipefail
[ "${ZSH_VERSION:-}" ] && setopt shwordsplit interactivecomments

: >log_DAJIN.txt
timestamp "Start analysis" log_DAJIN.txt
echo "DAJIN $*" >>log_DAJIN.txt
