#!/bin/bash

mamba create -n test-conda-DAJIN2 pytest --yes >/dev/null 2>&1
conda activate test-conda-DAJIN2
mamba install -n test-conda-DAJIN2 --file requirements.txt --yes || flag=1
python -m pytest -p no:warnings || flag=1
if [ "${flag:-0}" -eq 0 ]; then
    conda activate base
    mamba remove -n test-conda-DAJIN2 --all
else
    echo "TESTING WITH CONDA FAILED" && exit 1
fi

mamba create -n test-pypi-DAJIN2 >/dev/null 2>&1
conda activate test-pypi-DAJIN2
python -m pip install --upgrade pip
python -m pip install pytest
python -m pip install -r requirements.txt || flag=1
python -m pytest -p no:warnings || flag=1
if [ "${flag:-0}" -eq 0 ]; then
    conda activate base
    mamba remove -n test-pypi-DAJIN2 --all
else
    echo "TESTING WITH PIP FAILED" && exit 1
fi

echo "==================================="
echo "ALL TESTS ARE COMPLETED!"
echo "==================================="
