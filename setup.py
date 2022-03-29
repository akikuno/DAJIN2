from setuptools import setup

setup(
    name="DAJIN2",
    version="1.0.0",
    description="One-step genotyping tools for Nanopore amplicon sequencing",
    author="Akihiro Kuno",
    author_email="akuno@md.tsukuba.ac.jp",
    url="https://github.com/akikuno/DAJIN2",
    install_requires=["numpy", "pandas", "scikit-learn", "hdbscan", "plotnine", "mappy", "pysam"],
)
