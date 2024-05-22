#!/bin/bash

pip install -qU pip . || exit 1

PIP_VERSION=$(echo $(DAJIN2 --version) | awk '{print $NF}')
UTIL_VERSION=$(grep "DAJIN_VERSION" src/DAJIN2/utils/config.py | awk '{print $NF}' | tr -d '"')

if [ -z "$PIP_VERSION" ]; then
    echo "DAJIN2 is not installed"
    exit 1
fi

if [ -z "$UTIL_VERSION" ]; then
    echo "src/DAJIN2/utils/config.py is not found"
    exit 1
fi

if [ "$PIP_VERSION" != "$UTIL_VERSION" ]; then
    echo "Version mismatch: $PIP_VERSION != $UTIL_VERSION"
    exit 1
fi
