#!/bin/bash

sampath=".tmpDAJIN/sam/barcode30_target.sam"
extract_readid="e657b9f9-ffab-441c-96d8-2a15fbb95ca4"

cat $sampath |
    grep -e "^@" -e "$extract_readid" |
    samtools sort >tmp2.bam

samtools index tmp2.bam
