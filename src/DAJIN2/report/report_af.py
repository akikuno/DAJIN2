from __future__ import annotations
import pandas as pd
from collections import defaultdict
from plotnine import ggplot, aes, geom_bar, theme, theme_bw, element_blank, labs, scale_y_continuous


def call_allele_name(clust_sample: list[dict], cons_sequence: dict, dict_allele: dict) -> str:
    for c in clust_sample:
        _, ALLELE, SV, LABEL = c.values()
        allele_name = f"allele{LABEL}_{ALLELE}"
        key = f'{{"ALLELE": "{ALLELE}", "SV": {SV}, "LABEL": {LABEL}}}'
        if SV:
            allele_name += "_sv"
        elif cons_sequence[key] == dict_allele[ALLELE]:
            allele_name += "_intact"
        else:
            allele_name += "_mutated"
        c["NAME"] = allele_name
    return clust_sample


def all_allele(clust_sample: list[dict], sample_name: str) -> pd.DataFrame:
    df_clust_sample = pd.DataFrame(clust_sample)
    df_clust_sample["SAMPLE"] = sample_name
    colum = df_clust_sample.columns.to_list()
    colum = colum[-1:] + colum[:-1]
    return df_clust_sample[colum]


def summary_allele(clust_sample: list[dict], sample_name: str) -> pd.DataFrame:
    # #reads
    num_reads = defaultdict(int)
    for c in clust_sample:
        num_reads[c["LABEL"]] += 1
    # %reads
    per_reads = defaultdict(float)
    for LABEL, value in num_reads.items():
        per = value / len(clust_sample) * 100
        per_reads[LABEL] = float(f"{per:.2f}")
    # allele name
    allele_names = defaultdict(str)
    for c in clust_sample:
        if allele_names[c["LABEL"]]:
            continue
        allele_names[c["LABEL"]] = c["NAME"]
    allele_frequency = defaultdict(list)
    allele_frequency["sample"] = [sample_name] * len(num_reads)
    allele_frequency["allele name"] = [x for x in allele_names.values()]
    allele_frequency[r"#reads"] = [x for x in num_reads.values()]
    allele_frequency[r"%reads"] = [x for x in per_reads.values()]
    return pd.DataFrame(allele_frequency)


def plot(df_allele_frequency: pd.DataFrame):
    g = (
        ggplot(df_allele_frequency, aes(x="sample", y="#reads", fill="allele name"))
        + geom_bar(position="fill", stat="identity", colour="black")
        + scale_y_continuous(labels=lambda l: ["%d%%" % (v * 100) for v in l])
        + theme_bw()
        + theme(axis_title_x=element_blank())
        + labs(fill="Allele type", y="Percentage of reads")
    )
    return g
