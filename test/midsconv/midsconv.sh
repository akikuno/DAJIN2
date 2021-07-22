#!/bin/sh

. library/midsconv.sh

# >test/tmp_midsconv

find test/midsconv/input-*.sam |
  while read -r line; do
    output="$(echo "${line%.sam}".csv | sed "s/input/output/")"
    midsconv "$line" >test/tmp_midsconv
    diff test/tmp_midsconv "$output" ||
      echo "TEST FAILED: $output"
    rm test/tmp_midsconv
  done
