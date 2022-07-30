#!/bin/bash

cat .tmpDAJIN/midsv/barcode25_control.csv |
    sed "s|.*CSSPLIT||" |
    sed "s|QSCORE.*||" |
    head
