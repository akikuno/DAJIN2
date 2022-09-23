#!/bin/bash

mamba create -n test-conda-DAJIN2 pytest --yes
conda activate test-conda-DAJIN2
mamba install -n test-conda-DAJIN2 --file requirements.txt --yes || flag=1
python -m pytest -p no:warnings || flag=1
if [ "${flag:-0}" -eq 0 ]; then
    conda activate base
    mamba remove -n test-conda-DAJIN2 --all
fi

mamba create -n test-pypi-DAJIN2
conda activate test-pypi-DAJIN2
python -m pip install --upgrade pip
python -m pip install pytest
python -m pip install -r requirements.txt || flag=1
python -m pytest -p no:warnings || flag=1
if [ "${flag:-0}" -eq 0 ]; then
    conda activate base
    mamba remove -n test-pypi-DAJIN2 --all
fi
