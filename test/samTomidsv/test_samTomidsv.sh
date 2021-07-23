#!/bin/sh

. library/samTomidsv.sh

# >test/tmp_samTomidsv

ref=test/samTomidsv/test_ref.fa
ref_long=test/samTomidsv/test_ref_long.fa

find test/samTomidsv/*.fq |
  grep -v -e long -e inv |
  while read -r que; do
    minimap2 -ax map-ont "$ref" "$que" --cs=long 2>/dev/null >"${que%.fq}".sam
  done

find test/samTomidsv/*.fq |
  grep -e long -e inv |
  while read -r que; do
    minimap2 -ax map-ont "$ref_long" "$que" --cs=long 2>/dev/null >"${que%.fq}".sam
  done

find test/samTomidsv/input-*.sam |
  while read -r line; do
    output="$(echo "${line%.sam}".csv | sed "s/input/output/")"
    cat "$line" | samTomidsv >test/tmp_samTomidsv
    diff test/tmp_samTomidsv "$output" ||
      echo "TEST FAILED: $line"
    rm test/tmp_samTomidsv
  done
