#!/bin/sh

. library/midsconv.sh

# >test/tmp_midsconv

ref=test/midsconv/test_ref.fa
ref_long=test/midsconv/test_ref_long.fa

find test/midsconv/*.fq |
  grep -v -e long -e inv |
  while read -r que; do
    minimap2 -ax map-ont "$ref" "$que" --cs=long 2>/dev/null >"${que%.fq}".sam
  done

find test/midsconv/*.fq |
  grep -e long -e inv |
  while read -r que; do
    minimap2 -ax map-ont "$ref_long" "$que" --cs=long 2>/dev/null >"${que%.fq}".sam
  done

find test/midsconv/input-*.sam |
  while read -r line; do
    output="$(echo "${line%.sam}".csv | sed "s/input/output/")"
    cat "$line" | midsconv >test/tmp_midsconv
    diff test/tmp_midsconv "$output" ||
      echo "TEST FAILED: $line"
    rm test/tmp_midsconv
  done
