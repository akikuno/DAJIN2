#!/bin/sh

# find test/ -type f |
#   grep ".sh" |
#   xargs chmod +x

find test/ -type f |
  grep ".sh" |
  while read -r line; do
    if sh "$line"; then
      :
    else
      error at "$line"
    fi
  done
