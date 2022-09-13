from __future__ import annotations
from src.DAJIN2.report import report_af
import pandas as pd
from importlib import reload
from pathlib import Path

reload(report_af)


def test_plot():
    df = pd.read_csv("tests/plot_alleles/test_input.csv")
    g_test = report_af.plot(df)
    g_test.save(filename="tests/plot_alleles/test.eps")
    test = Path("tests/plot_alleles/test.eps").read_text().split("\n")
    test = [t for t in test if not t.startswith(r"%")]
    answer = Path("tests/plot_alleles/answer.eps").read_text().split("\n")
    answer = [a for a in answer if not a.startswith(r"%")]
    assert test == answer
