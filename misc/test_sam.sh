#!/bin/bash

sampath=".tmpDAJIN/sam/barcode31_albino.sam"

cat $sampath |
    grep -e "^@" -e bdf684fb-2ae6-4a0e-9fee-73d30778c457 |
    samtools sort >tmp.bam

samtools index tmp.bam
