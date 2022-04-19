[![licence](https://img.shields.io/badge/License-MIT-blue.svg?style=flat-square)](https://choosealicense.com/licenses/mit/)
<!-- [![PyPI version](https://img.shields.io/badge/Install%20with-PyPI-brightgreen.svg?style=flat-square)](https://pypi.org/project/DAJIN2/)
[![install with bioconda](https://img.shields.io/badge/Install%20with-Bioconda-brightgreen.svg?style=flat-square)](https://anaconda.org/bioconda/DAJIN2) -->

# DAJIN2

## Description

## Installation

You can install `DAJIN2` using bioconda or pip:

```bash
conda install -c bioconda DAJIN2
# conda create -y -n DAJIN2
# conda install -y -n DAJIN2 cstag numpy scikit-learn hdbscan plotnine mappy pysam
# conda activate DAJIN2
```

```bash
pip install DAJIN2
```

<!-- > :warning: install `minimap2` and `samtools` when you use pip. -->

## Usage

```bash
DAJIN [options] -a/allele <alleles.fa> -c <control.fasta> -s <sample.fasta>
```

## Options

```bash
-o/--output [STR]: output directory <default: DAJIN2_results>

-g/--genome [STR]: UCSC Genome assembly ID. (e.g. hg38, mm10)

-t/--threads [INT]: number of threads to use <default: 1>
```

## Examples

```bash
:
```