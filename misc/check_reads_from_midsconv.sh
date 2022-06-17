#!/bin/bash

sample=barcode25
readid_overlapped=f44a2a8c-5a9a-4b94-86da-7f19529cc8df #<- overlapped reads in barcode25, stx2
readid_largedeletion=cd66bf2a-672b-42b1-8139-023f201476fc

minimap2 -ax map-ont -t 20 tests/data/mappy/stx2-wt.fa examples/del-stx2/barcode25.fq.gz |
    grep -e "^@" -e "$readid_overlapped" -e "$readid_largedeletion" |
    cat >overlapped.sam

# minimap2 -ax map-ont -t 20 .tmpDAJIN/fasta/control.fasta examples/del-stx2/barcode25.fq.gz --cs=long> tmp_largedel.sam
# samtools sort tmp_largedel.sam > tmp_largedel.bam
# samtools index tmp_largedel.bam
# cat tmp_largedel.sam | grep 023f201476fc | cut -f 1-6 #<- large deletion

# minimap2 -ax map-ont -t 20 .tmpDAJIN/fasta/inversion.fasta examples/del-stx2/barcode30.fq.gz > tmp.sam
minimap2 -ax map-ont -t 20 .tmpDAJIN/fasta/control.fasta examples/del-stx2/"$sample".fq.gz >tmp.sam
cat tmp.sam | grep "$readid" | cut -f 1-10

find .tmpDAJIN/midsconv/"$sample"* |
    while read -r line; do
        echo $line
        grep "$readid" "$line"
    done

cat .tmpDAJIN/sam/barcode25_target.sam |
    grep -e "^@" -e "$readid" |
    samtools sort >tmp.bam
samtools index tmp.bam

cat .tmpDAJIN/sam/barcode25_control.sam |
    grep -e "^@" -e "$readid" |
    samtools sort >tmp2.bam
samtools index tmp2.bam

# Large deletion

minimap2 -ax map-ont -t 20 .tmpDAJIN/fasta/control.fasta examples/del-stx2/barcode25.fq.gz >tmp.sam
cat tmp.sam | cut -f 1 | sort | uniq -c | awk '$1 == 2 {print $2}' >tmp_largedel
grep -e "^@" -f tmp_largedel tmp.sam >tmp_largedel.sam
samtools sort tmp_largedel.sam >tmp_largedel.bam
samtools index tmp_largedel.bam
