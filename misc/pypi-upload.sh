#!/bin/bash

pip uninstall DAJIN2 --yes
rm -rf dist build src/DAJIN2.egg-info
python3 -m pip install --upgrade pip setuptools wheel twine build
python3 setup.py sdist bdist_wheel

cat <<TESTPYPI
python3 -m twine upload --repository testpypi dist/*
TESTPYPI

python3 -m twine upload dist/*

# https://pypi.org/project/DAJIN2/
cat <<PIPINSTALL
pip uninstall DAJIN2 --yes && pip install DAJIN2
PIPINSTALL
