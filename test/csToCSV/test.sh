#!/bin/sh

minimap2 --cs=long -ax map-ont test/csToCSV/ref.fa test/csToCSV/que.fq >test/csToCSV/que.sam
minimap2 --cs=long -ax map-ont test/csToCSV/ref_long.fa test/csToCSV/que_long.fq >test/csToCSV/que_long.sam

head test/csToCSV/que.sam
head test/csToCSV/que_long.sam

. library/consensus/csToCSV.sh
cat test/csToCSV/que.sam |
  csToCSV
