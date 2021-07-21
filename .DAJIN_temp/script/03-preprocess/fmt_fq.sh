#!/bin/sh

fmt_fq() (
  {
    input="$1"
    if file "$input" | grep -q "compressed"; then
      CMD_CAT="gzip -dc"
    else
      CMD_CAT="cat"
    fi

    $CMD_CAT "$input" |
      tr -d "\r"
  } 2>.DAJIN_temp/error
  if [ -s .DAJIN_temp/error ]; then
    error_exit "$input" is not fastq file
  fi
)
