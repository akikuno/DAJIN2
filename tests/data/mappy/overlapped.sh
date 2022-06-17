#!/bin/bash

readid_overlapped=f44a2a8c-5a9a-4b94-86da-7f19529cc8df
readid_largedeletion=cd66bf2a-672b-42b1-8139-023f201476fc

minimap2 -ax map-ont -t 20 tests/data/mappy/stx2-wt.fa examples/del-stx2/barcode25.fq.gz |
    grep -e "^@" -e "$readid_overlapped" -e "$readid_largedeletion" |
    cat >tests/data/mappy/overlapped.sam

samtools view -h examples/flox-cables2/AyabeTask2/mm39/barcode47.bam |
    grep -e "^@SQ" -e ba5a5de6-c7e1-44b3-b7e1-86873435c268 >tests/data/mappy/inversion.sam

samtools view -h examples/flox-cables2/AyabeTask2/mm39/barcode47.bam >tmp.sam

cat tmp_remove_overlapped_reads.sam | cut -f 1 | sort -u >tmp1
cat tmp.sam | cut -f 1 | sort -u >tmp2

join -v 1 tmp2 tmp1 >tmp3

output="tmp_inversion"
samtools view -h examples/flox-cables2/AyabeTask2/mm39/barcode47.bam |
    grep -e "@SQ" -e 9c21963f917c >"$output".sam
samtools sort "$output".sam >"$output".bam
samtools index "$output".bam
