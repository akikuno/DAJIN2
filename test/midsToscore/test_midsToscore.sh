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

cat <<EOF |
aaa,M,2M,M,D,M
bbb,M,2M,D,D,M
ccc,M,1M,M,M,S
EOF
  midsToscore

cat <<EOF |
aaa,M,2M,M,D,M,M,M,M,M
bbb,M,2M,D,D,M,D,D,D,M
ccc,M,1M,M,D,S,S,M,M,M
EOF
  midsToscore

cat <<EOF |
aaa,D,D,M,D,D,M
bbb,D,D,M,D,D,M
EOF
  midsToscore
