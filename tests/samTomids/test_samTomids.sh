#!/bin/sh

. library/samTomids.sh

type minimap2

ref=test/samTomids/misc/test_ref.fa
ref_long=test/samTomids/misc/test_ref_long.fa

find test/samTomids/misc/*.fq |
  grep -v -e long -e inv |
  while read -r que; do
    minimap2 -ax map-ont "$ref" "$que" --cs=long >"${que%.fq}".sam
    mv "${que%.fq}".sam test/samTomids/in/
  done

find test/samTomids/misc/*.fq |
  grep -e long -e inv |
  while read -r que; do
    minimap2 -ax map-ont "$ref_long" "$que" --cs=long >"${que%.fq}".sam
    mv "${que%.fq}".sam test/samTomids/in/
  done

find test/samTomids/in/test*sam |
  while read -r line; do
    samTomids "$line" >"${line%.sam}.csv"
    mv "${line%.sam}.csv" test/samTomids/out/
  done

. library/samTomids.sh
cat test/samTomids/in/test_del.sam | samTomids
samTomids test/samTomids/in/test_del.sam
samTomids test/samTomids/in/test_padding.sam
samTomids test/samTomids/in/test_inv.sam

# find test/samTomids/input-*.sam |
#   while read -r line; do
#     output="$(echo "${line%.sam}".csv | sed "s/input/output/")"
#     cat "$line" | samTomids >test/tmp_samTomids
#     diff test/tmp_samTomids "$output" ||
#       echo "TEST FAILED: $line"
#     rm test/tmp_samTomids
#   done
