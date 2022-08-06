#!/bin/bash

suffix=$(date | tr " :" -)

sampath=".tmpDAJIN/sam/barcode31_flox.sam"
# sampath=".tmpDAJIN/sam/barcode42_flox.sam"
# extract_readid="7cdb4acdbb1b"
# extract_readid=ffbcdc27-95dd-4635-94bf-45c38c6771b1

# cat $sampath |
#     grep -e "^@" -e "${extract_readid-:*}" |
#     samtools sort >tmp_"$suffix".bam
# samtools index tmp_"$suffix".bam

filepath="tmp_qnames_1750_deletion.csv"
# filepath="tmp_qnames_2463.csv"
cat $sampath |
    grep -f "$filepath" |
    samtools sort >tmp_"$suffix".bam
samtools index tmp_"$suffix".bam
