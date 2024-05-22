#!/bin/bash

PIP_VERSION=$(echo $(DAJIN2 --version) | awk '{print $NF}')
UTIL_VERSION=$(grep "DAJIN_VERSION" DAJIN2/src/DAJIN2/utils/config.py | awk '{print $NF}' | tr -d '"')

if [ "$PIP_VERSION" != "$UTIL_VERSION" ]; then
    echo "Version mismatch: $PIP_VERSION != $UTIL_VERSION"
    exit 1
fi
