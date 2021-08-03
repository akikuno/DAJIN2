#!/bin/sh

. library/maskFastq.sh

cat test/maskFastq/test.fq
cat test/maskFastq/test.fq | maskFastq
