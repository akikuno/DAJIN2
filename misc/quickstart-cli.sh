#!/bin/bash

pip install -e . && DAJIN2 batch -f examples/flox-cables2/AyabeTask1/batch.csv --debug

pip install -e . && DAJIN2 gui --debug

pip install -e . && DAJIN2 view -n Ayabe-Task1
