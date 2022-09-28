#!/bin/bash

echo "Creating ENV"
mamba create -n test-conda-DAJIN2 pytest --yes >/dev/null 2>&1
mamba create -n test-pypi-DAJIN2 --yes >/dev/null 2>&1

cat <<EOF
===================================
TESTING WITH CONDA
===================================
EOF

conda activate test-conda-DAJIN2 >/dev/null 2>&1
mamba install -n test-conda-DAJIN2 --file requirements.txt --yes >/dev/null 2>&1 || flag=1
python -m pytest -p no:warnings || flag=1
if [ "${flag:-0}" -eq 0 ]; then
    conda activate base >/dev/null 2>&1
    mamba remove -n test-conda-DAJIN2 --all >/dev/null 2>&1
else
    echo "TESTING WITH CONDA FAILED" && exit 1
fi

cat <<EOF
===================================
TESTING WITH PYPI
===================================
EOF
conda activate test-pypi-DAJIN2 >/dev/null 2>&1
python -m pip install --upgrade pip >/dev/null 2>&1
python -m pip install pytest >/dev/null 2>&1
python -m pip install -r requirements.txt >/dev/null 2>&1 || flag=1
python -m pytest -p no:warnings || flag=1
if [ "${flag:-0}" -eq 0 ]; then
    conda activate base >/dev/null 2>&1
    mamba remove -n test-pypi-DAJIN2 --all >/dev/null 2>&1
else
    echo "TESTING WITH PIP FAILED" && exit 1
fi

cat <<EOF
===================================
ALL TESTS ARE COMPLETED!
===================================
EOF
