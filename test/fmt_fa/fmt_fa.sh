#!/bin/sh

. library/fmt_fa.sh

fmt_fa test/fmt_fa/in1.fa >tmp1
diff tmp1 test/fmt_fa/out1
fmt_fa test/fmt_fa/in2.fa >tmp2
diff tmp2 test/fmt_fa/out2
