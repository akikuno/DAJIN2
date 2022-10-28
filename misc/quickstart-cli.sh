#!/bin/bash

pip uninstall DAJIN2 -y && pip install -e .

rm -rf DAJINResults/test-single DAJINResults/.tempdir/test-single
pip install -e . && DAJIN2 \
    --name test-single \
    --sample misc/data/tyr_albino_01%.fq.gz \
    --control misc/data/tyr_control.fq.gz \
    --allele misc/data/tyr_control.fasta
DAJIN2 view -n test-single

rm -rf DAJINResults/test-single DAJINResults/.tempdir/test-single
pip install -e . && DAJIN2 \
    --name test-single \
    --sample examples/flox-cables2/AyabeTask1/barcode32.fq.gz \
    --control examples/flox-cables2/AyabeTask1/barcode42.fq.gz \
    --allele examples/flox-cables2/AyabeTask1/design_cables2.fa

rm -rf DAJINResults/Ayabe-Task1 DAJINResults/.tempdir/Ayabe-Task1
pip install -e . && DAJIN2 batch -f examples/flox-cables2/AyabeTask1/batch.csv --debug

pip install -e . && DAJIN2 gui --debug

pip install -e . && DAJIN2 view -n Ayabe-Task1
