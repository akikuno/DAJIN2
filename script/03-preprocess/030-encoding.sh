#!/bin/sh

#----------------------------------------------------------
echo "$(date +'%Y-%m-%d %H:%M:%S') MIDS conversion" >>log_DAJIN.txt
#----------------------------------------------------------

mkdir -p .DAJIN_temp/mids
. .DAJIN_temp/library/midsconv.sh

cmd='. .DAJIN_temp/library/midsconv.sh && midsconv'

find .DAJIN_temp/sam/*.sam |
  grep -e "$control_name" -e "$sample_name" |
  awk -v th="$threads" -v cmd="$cmd" '{
    input=$0; output=$0
    sub("sam$", "csv", output)
    sub(".DAJIN_temp/sam/", ".DAJIN_temp/mids/", output)
    print "[ -f ",output, "] || {",cmd, input, ">", output, ";} &"
    if (NR%th==0) print "wait"
    } END {
    print "wait"
  }' |
  sh

# Error handling ------------------------------------------
find .DAJIN_temp/mids/ -type f |
  while read -r line; do
    if ! [ -s "${line}" ]; then
      error_exit "${line} is empty at MIDS conversion"
      break
    fi
  done || exit 1
