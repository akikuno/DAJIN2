#!/bin/sh

cat <<EOF >tmp_s.csv
a,D,D,D
b,D,D,D
c,D,D,M
d,D,S,M
EOF

cat <<EOF >tmp_c.csv
a,D,D,D
b,D,D,D
c,D,D,M
d,D,S,M
e,D,M,M
EOF

set tmp_s.csv tmp_c.csv
set .DAJIN_temp/midsmask/barcode31_control.csv .DAJIN_temp/midsmask/barcode32_control.csv
cp -f $1 tmp_s.csv
cp -f $2 tmp_c.csv
fn="$(find library -name "mids-subtract.R")"
time Rscript --vanilla --slave "$fn" "$1" "$2" >tmp
