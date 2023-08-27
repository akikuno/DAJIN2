[![License](https://img.shields.io/badge/License-MIT-9cf.svg?style=flat-square)](https://choosealicense.com/licenses/mit/)
[![Test](https://img.shields.io/github/actions/workflow/status/akikuno/dajin2/pytest.yml?branch=main&label=Test&color=brightgreen&style=flat-square)](https://github.com/akikuno/dajin2/actions)
[![Python](https://img.shields.io/pypi/pyversions/DAJIN2.svg?label=Python&color=blue&style=flat-square)](https://pypi.org/project/DAJIN2/)
[![PyPI](https://img.shields.io/pypi/v/DAJIN2.svg?label=PyPI&color=orange&style=flat-square)](https://pypi.org/project/DAJIN2/)
[![Bioconda](https://img.shields.io/conda/v/bioconda/dajin2?label=Bioconda&color=orange&style=flat-square)](https://anaconda.org/bioconda/dajin2)


<p align="center">
<img src="https://user-images.githubusercontent.com/15861316/261833016-7f356960-88cf-4574-87e2-36162b174340.png" width="90%">
</p>

DAJIN2 is a genotyping software designed for organisms that have undergone genome editing, utilizing nanopore sequencing technology.  

The name DAJIN is inspired by the term 一網**打尽** (Ichimou **DAJIN** or Yīwǎng **Dǎjìn**), which signifies capturing everything in a single net.  

## Disclaimer

DAJIN2 is still in the development phase.  
Basic tests covering point mutations, deletions, and insertion designs have been conducted.  
If you encounter any bugs or issues, please report them via [Issues](https://github.com/akikuno/DAJIN2/issues).  



## Installation

From [PyPI](https://pypi.org/project/DAJIN2/):

```bash
pip install DAJIN2
```

From [Bioconda](https://anaconda.org/bioconda/DAJIN2):

```bash
conda install -c bioconda DAJIN2
```


## Usage

### Single Sample Analysis

DAJIN2 allows for the analysis of single samples (one sample vs one control).

```bash
DAJIN2 <-s|--sample> <-c|--control> <-a|--allele> <-n|--name> [-g|--genome] [-t|--threads] [-h|--help] [-v|--version]

options:
  -s, --sample              Path to a sample FASTQ file
  -c, --control             Path to a control FASTQ file
  -a, --allele              Path to a FASTA file
  -n, --name                Output directory name
  -g, --genome (Optional)   Reference genome ID (e.g hg38, mm39) [default: '']
  -t, --threads (Optional)  Number of threads [default: 1]
  -h, --help                show this help message and exit
  -v, --version             show the version number and exit
```

#### Example

```bash
# Donwload the example dataset
wget https://github.com/akikuno/DAJIN2/raw/main/examples/example-single.tar.gz
tar -xf example-single.tar.gz

# Run DAJIN2
DAJIN2 \
    --name stx2-deletion \
    --sample example-single/sample.fq.gz \
    --control example-single/control.fq.gz \
    --allele example-single/design.fa \
    --genome mm39 \
    --threads 10

# 2023-06-04 11:30:03: example-single/control.fq.gz is now processing...
# 2023-06-04 11:30:06: Preprocess example-single/control.fq.gz...
# 2023-06-04 11:30:06: Mapping example-single/control.fq.gz...
# 2023-06-04 11:30:21: Call MIDSV example-single/control.fq.gz...
# 2023-06-04 11:30:31: 🍵 example-single/control.fq.gz is finished!
# 2023-06-04 11:30:31: example-single/sample.fq.gz is now processing...
# 2023-06-04 11:30:35: Preprocess example-single/sample.fq.gz...
# 2023-06-04 11:34:13: Classify example-single/sample.fq.gz...
# 2023-06-04 11:34:18: Clustering example-single/sample.fq.gz...
# 2023-06-04 11:35:01: Consensus calling example-single/sample.fq.gz...
# 2023-06-04 11:35:08: 🍵 example-single/sample.fq.gz is finished!
# 🎉 Finished! Open DAJINResults/stx2-deletion to see the report.
```

### Batch Processing

DAJIN2 can also handle multiple FASTQ files using the `batch` subcommand.


```bash
DAJIN2 batch <-f|--file> [-t|--threads] [-h]

options:
  -f, --file                Path to a CSV or Excel file
  -t, --threads (Optional)  Number of threads [default: 1]
  -h, --help                Show this help message and exit
```

#### Example

```bash
# Donwload the example dataset
wget https://github.com/akikuno/DAJIN2/raw/main/examples/example-batch.tar.gz
tar -xf example-batch.tar.gz

# Run DAJIN2
DAJIN2 batch --file example-batch/batch.csv --threads 3

# 2023-07-31 17:01:10: example-batch/tyr_control.fq.gz is now processing...
# 2023-07-31 17:01:16: Preprocess example-batch/tyr_control.fq.gz...
# 2023-07-31 17:01:48: Output BAM files of example-batch/tyr_control.fq.gz...
# 2023-07-31 17:01:52: 🍵 example-batch/tyr_control.fq.gz is finished!
# 2023-07-31 17:01:52: example-batch/tyr_c230gt_50%.fq.gz is now processing...
# 2023-07-31 17:01:52: example-batch/tyr_c230gt_10%.fq.gz is now processing...
# 2023-07-31 17:01:52: example-batch/tyr_c230gt_01%.fq.gz is now processing...
# 2023-07-31 17:01:55: Preprocess example-batch/tyr_c230gt_01%.fq.gz...
# 2023-07-31 17:01:55: Preprocess example-batch/tyr_c230gt_50%.fq.gz...
# 2023-07-31 17:01:55: Preprocess example-batch/tyr_c230gt_10%.fq.gz...
# 2023-07-31 17:02:17: Classify example-batch/tyr_c230gt_50%.fq.gz...
# 2023-07-31 17:02:19: Clustering example-batch/tyr_c230gt_50%.fq.gz...
# 2023-07-31 17:02:34: Classify example-batch/tyr_c230gt_01%.fq.gz...
# 2023-07-31 17:02:35: Classify example-batch/tyr_c230gt_10%.fq.gz...
# 2023-07-31 17:02:39: Clustering example-batch/tyr_c230gt_01%.fq.gz...
# 2023-07-31 17:02:39: Clustering example-batch/tyr_c230gt_10%.fq.gz...
# 2023-07-31 17:02:53: Consensus calling of example-batch/tyr_c230gt_50%.fq.gz...
# 2023-07-31 17:02:59: Output reports of example-batch/tyr_c230gt_50%.fq.gz...
# 2023-07-31 17:03:04: 🍵 example-batch/tyr_c230gt_50%.fq.gz is finished!
# 2023-07-31 17:03:39: Consensus calling of example-batch/tyr_c230gt_01%.fq.gz...
# 2023-07-31 17:03:51: Output reports of example-batch/tyr_c230gt_01%.fq.gz...
# 2023-07-31 17:04:03: 🍵 example-batch/tyr_c230gt_01%.fq.gz is finished!
# 2023-07-31 17:04:08: Consensus calling of example-batch/tyr_c230gt_10%.fq.gz...
# 2023-07-31 17:04:16: Output reports of example-batch/tyr_c230gt_10%.fq.gz...
# 2023-07-31 17:04:24: 🍵 example-batch/tyr_c230gt_10%.fq.gz is finished!
# 🎉 Finished! Open DAJINResults/tyr-substitution to see the report.
```

## References

For more information, please refer to the following publication:

[Kuno A, et al. (2022) DAJIN enables multiplex genotyping to simultaneously validate intended and unintended target genome editing outcomes. *PLoS Biology* 20(1): e3001507.](https://doi.org/10.1371/journal.pbio.3001507)
