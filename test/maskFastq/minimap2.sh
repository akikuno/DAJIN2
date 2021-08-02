#!/bin/sh


minimap2 -ax map-ont ref.fa test_mask.fq --cs=long > tmp.sam

cat tmp.sam