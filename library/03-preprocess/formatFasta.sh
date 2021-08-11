#!/bin/sh

formatFasta() {
  cat "$1" |
    tr -d "\r" |
    paste - - |
    awk 'NF==2 {$2=toupper($2)}1'
}
