from __future__ import annotations
from src.DAJIN2.report import report_af
import pandas as pd
from importlib import reload
from pathlib import Path
from collections import defaultdict
import tempfile

reload(report_af)


def test_all_allele():
    x = [
        {"QNAME": "nobita", "ALLELE": "albino", "SV": False, "LABEL": 1},
        {"QNAME": "shizuka", "ALLELE": "albino", "SV": False, "LABEL": 2},
        {"QNAME": "suneo", "ALLELE": "control", "SV": False, "LABEL": 3},
        {"QNAME": "gian", "ALLELE": "control", "SV": True, "LABEL": 4},
        {"QNAME": "dora", "ALLELE": "albino", "SV": True, "LABEL": 5},
    ]
    df_clust_sample = pd.DataFrame(x)
    df_clust_sample["SAMPLE"] = "test"
    colum = df_clust_sample.columns.to_list()
    colum = colum[-1:] + colum[:-1]
    test = df_clust_sample[colum]
    answer = pd.read_csv("tests/data/report_af_all_allele/answer.csv")
    assert test.to_json() == answer.to_json()


def test_summary_allele():
    # prepare inputs
    clust_sample = Path("tests/data/report_af_summary_allele/clust_sample.txt").read_text()
    clust_sample = eval(clust_sample)
    cons_sequence = Path("tests/data/report_af_summary_allele/cons_sequence.txt").read_text()
    cons_sequence = eval(cons_sequence)
    dict_allele = Path("tests/data/report_af_summary_allele/dict_allele.txt").read_text()
    dict_allele = eval(dict_allele)
    sample_name = "barcode31"
    test = report_af.summary_allele(clust_sample, sample_name, cons_sequence, dict_allele)
    answer = pd.read_csv("tests/data/report_af_summary_allele/answer.csv")
    assert test.to_json() == answer.to_json()


# def test_plot():
#     df = pd.read_csv("tests/data/report_af_plot/test_input.csv")
#     g_test = report_af.plot(df)
#     filename = tempfile.NamedTemporaryFile().name + ".eps"
#     g_test.save(filename=filename)
#     test = Path(filename).read_text().split("\n")
#     test = [t for t in test if not t.startswith(r"%")]
#     answer = Path("tests/data/report_af_plot/answer.eps").read_text().split("\n")
#     answer = [a for a in answer if not a.startswith(r"%")]
#     assert test == answer
