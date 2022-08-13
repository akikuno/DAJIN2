#!/bin/bash

suffix=$(date | tr " :" -)

# sampath=".tmpDAJIN/sam/barcode31_flox.sam"
# sampath=".tmpDAJIN/sam/barcode42_flox.sam"
sampath=".tmpDAJIN/sam/barcode31_control.sam"

# extract_readid="6f536ce8-32d6-41d5-8b09-66bef425b9d8"
# cat $sampath |
#     grep -e "^@" -e "${extract_readid-:*}" |
#     samtools sort >tmp_"$suffix".bam
# samtools index tmp_"$suffix".bam

# cat <<EOF >/tmp/tmp_reads.txt
# ^@
# *
# EOF

# filepath=/tmp/tmp_reads.txt
# cat $sampath |
#     grep -f "$filepath" |
#     samtools sort >tmp_barcode31_control_"$suffix".bam
# samtools index tmp_barcode31_control_"$suffix".bam

filepath="tmp_0"
cat $sampath |
    grep -f "$filepath" |
    samtools sort >"$filepath"_"$suffix".bam
samtools index "$filepath"_"$suffix".bam
