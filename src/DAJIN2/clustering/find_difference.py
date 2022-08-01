from pathlib import Path
import scipy.stats as st
import midsv

ALLELE = "target"

sample_cssplit = Path(".tmpDAJIN", "midsv", f"{sample_name}_{ALLELE}.csv")
sample_cssplit = midsv.read_jsonl(sample_cssplit)
sample_cssplit = [cs["CSSPLIT"] for cs in sample_cssplit]

control_cssplit = Path(".tmpDAJIN", "midsv", f"{control_name}_{ALLELE}.csv")
control_cssplit = midsv.read_jsonl(control_cssplit)
control_cssplit = [cs["CSSPLIT"] for cs in control_cssplit]

length = len(control_cssplit[0].split(","))
sample_table = [[0, 0] for _ in range(length)]
control_table = [[0, 0] for _ in range(length)]

for sample_cs, control_cs in zip(sample_cssplit, control_cssplit):
    sample_cs = sample_cs.split(",")
    control_cs = control_cs.split(",")
    for i, (s, c) in enumerate(zip(sample_cs, control_cs)):
        if s.startswith("="):
            sample_table[i][0] += 1
        else:
            sample_table[i][1] += 1
        if c.startswith("="):
            control_table[i][0] += 1
        else:
            control_table[i][1] += 1

s = sample_table[0]
c = control_table[0]

s = sample_table[400]
c = control_table[400]
st.chi2_contingency([s, c])[:2]

table_test = []
for i, (s, c) in enumerate(zip(sample_table, control_table)):
    odds = (s[1] * c[0]) / (s[0] * c[1])
    pval = st.chi2_contingency([s, c])[1]
    if odds > 1 and pval < 0.01:
        table_test.append((i, odds, pval))

table_test
len(table_test)

