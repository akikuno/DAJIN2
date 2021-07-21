#!/bin/sh

set -eu
umask 0022
export LC_ALL=C

[ "${BASH_VERSION:-}" ] && set -o posix -o pipefail
[ "${ZSH_VERSION:-}" ] && setopt shwordsplit interactivecomments

timestamp "DAJIN $*" >log_DAJIN.txt
