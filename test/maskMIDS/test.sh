#!/bin/sh

. library/samTomisa
minimap2 -ax map-ont --cs=long test/maskMIDS/ref.fa test/maskMIDS/que.fq >tmp.sam

cat tmp.sam |

