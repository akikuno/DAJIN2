#!/bin/bash

rm -rf dist build DAJIN2.egg-info
pip uninstall DAJIN2 --yes
python3 -m pip install --upgrade setuptools wheel twine build
python3 setup.py sdist bdist_wheel
python3 -m twine upload dist/*

# https://pypi.org/project/DAJIN2/
cat <<MEMO
pip install DAJIN2==0.1.3
MEMO
