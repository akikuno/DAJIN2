import csv
import warnings
import pandas as pd

from plotnine import ggplot, aes, geom_bar, theme, theme_bw, element_blank, labs

warnings.simplefilter("ignore")


file = "tests/plot_alleles/test_input.csv"
with open(file) as f:
    reader = csv.reader(f)
    input = [row for row in reader]

sample = [s for s, a, r in input]
allele = [a for s, a, r in input]
readnum = [float(r) for s, a, r in input]

df = pd.DataFrame({"sample": sample, "allele": allele, "readnum": readnum})

g = (
    ggplot(df, aes(x="sample", y="readnum", fill="allele"))
    + geom_bar(position="fill", stat="identity", colour="black")
    + theme_bw()
    + theme(axis_title_x=element_blank())
    + labs(fill="Allele type", y="Percentage of Alleles")
)

g.save(filename="tests/plot_alleles/test_output.png", dpi=350)
