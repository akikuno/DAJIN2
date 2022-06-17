#!/bin/bash

readid_overlapped=f44a2a8c-5a9a-4b94-86da-7f19529cc8df
readid_largedeletion=cd66bf2a-672b-42b1-8139-023f201476fc

minimap2 -ax map-ont -t 20 tests/data/mappy/stx2-wt.fa examples/del-stx2/barcode25.fq.gz |
    grep -e "^@" -e "$readid_overlapped" -e "$readid_largedeletion" |
    cat >tests/data/mappy/overlapped.sam
