#!/bin/bash

ref="tests/data/mappy/ref.fa"
que="tests/data/mappy/query.fq"

minimap2 -ax map-ont "$ref" "$que" >tests/data/mappy/align.sam
minimap2 -ax map-ont "$ref" "$que" --cs >tests/data/mappy/align_cs.sam
minimap2 -ax map-ont "$ref" "$que" --cs=long >tests/data/mappy/align_cslong.sam

ref=".tmpDAJIN/fasta/control.fasta"
que="examples/pm-tyr/barcode31.fq.gz"
threads=14

minimap2 -t "$threads" -ax map-ont "$ref" "$que" >tmp_minimap2.sam
samtools sort -@ "$threads" tmp_minimap2.sam >tmp_minimap2.bam
samtools index -@ "$threads" tmp_minimap2.bam

id="accf9e5ceb7a"
cat tmp.sam | grep "$id" | cut -f 1-10
cat tmp_minimap2.sam | grep "$id" | cut -f 1-10
