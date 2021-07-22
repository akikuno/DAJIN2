#!/bin/sh

. library/fmt_fa.sh

fmt_fa test/fmt_fa/in1.fa >/tmp/tmp1
diff /tmp/tmp1 test/fmt_fa/out1
fmt_fa test/fmt_fa/in2.fa >/tmp/tmp2
diff /tmp/tmp2 test/fmt_fa/out2

rm /tmp/tmp*
