#!/bin/bash

sampath=".tmpDAJIN/sam/barcode31_albino.sam"

cat $sampath |
    grep -e "^@" -e 001179a7-d376-42d5-9671-dde7bff87c93 |
    samtools sort >tmp.bam

samtools index tmp.bam
