#!/bin/sh

python src/DAJIN2/main.py \
  -d -t 20 \
  -s tests/data/query.fq \
  -c tests/data/query.fq \
  -a tests/data/ref.fa
