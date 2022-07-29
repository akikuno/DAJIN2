#!/bin/bash

suffix=$(date | tr " " -)

sampath=".tmpDAJIN/sam/barcode25_control.sam"
extract_readid="bfcf883e9e86"
# extract_readid=ffbcdc27-95dd-4635-94bf-45c38c6771b1

cat $sampath |
    grep -e "^@" -e "$extract_readid" |
    samtools sort >tmp_"$suffix".bam

samtools index tmp_"$suffix".bam
