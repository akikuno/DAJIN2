#!/bin/sh

. library/midsToscore.sh

echo "aaa,M,M,M,M" |
  midsToscore

echo "aaa,S,S,S,S,S" |
  midsToscore

echo "aaa,M,100S,D,M" |
  midsToscore

echo "aaa,D,100S,D,D" |
  midsToscore
