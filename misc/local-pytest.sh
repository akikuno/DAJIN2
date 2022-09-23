#!/bin/bash

mamba create -n test-DAJIN2 pytest --yes
mamba install -n test-DAJIN2 --file requirements.txt --yes
conda activate test-DAJIN2
if pytest; then
    conda activate base
    mamba remove -n test-DAJIN2 --all
fi
