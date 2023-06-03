#!/bin/bash

python -m pip install --upgrade pip pyOpenSSL

python -m pip uninstall DAJIN2 --yes
rm -rf dist build src/DAJIN2.egg-info
python -m pip install --upgrade pip setuptools wheel twine build
python setup.py sdist bdist_wheel

cat <<TESTPYPI
python -m twine upload --repository testpypi dist/*
TESTPYPI

python -m twine upload dist/*

# https://pypi.org/project/DAJIN2/
cat <<PIPINSTALL
pip install -U DAJIN2
PIPINSTALL
