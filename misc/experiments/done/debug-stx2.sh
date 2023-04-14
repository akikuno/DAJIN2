#!/bin/bash

cat examples/del-stx2/design_stx2.fa | grep "control" -A1 >tmpstx.fa
zcat examples/del-stx2/barcode25.fq.gz | grep 3f7aff2243b8 -A3 >tmp3f7aff2243b8.fq
minimap2 -ax map-ont tmpstx.fa tmp3f7aff2243b8.fq --cs=long >tmp3f7aff2243b8.sam

cat DAJINResults/.tempdir/stx2-deletion/sam/barcode25_control.sam | grep 3f7aff2243b8 >tmp3f7aff2243b8_DAJIN2.sam
