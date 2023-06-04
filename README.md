[![License](https://img.shields.io/badge/License-MIT-9cf.svg?style=flat-square)](https://choosealicense.com/licenses/mit/)
[![PyPI](https://img.shields.io/pypi/v/DAJIN2.svg?label=PyPI&color=orange&style=flat-square)](https://pypi.org/project/DAJIN2/)

⚠️ DAJIN2 is currently under development ⚠️

Expected to be available the stable version in August 2023 🤞

## Installation (alpha-version)

```bash
pip install DAJIN2
```

## Usage

### Single-mode

```bash
DAJIN2 [-h] [-s SAMPLE] [-c CONTROL] [-a ALLELE] [-n NAME] [-g GENOME] [-t THREADS] [-v]

options:
  -h, --help            show this help message and exit
  -s SAMPLE, --sample SAMPLE
                        Full path to a sample FASTQ file
  -c CONTROL, --control CONTROL
                        Full path to a control FASTQ file
  -a ALLELE, --allele ALLELE
                        Full path to a FASTA file
  -n NAME, --name NAME  Output directory name
  -g GENOME, --genome GENOME
                        Reference genome ID (e.g hg38, mm10) [default: '']
  -t THREADS, --threads THREADS
                        Number of threads [default: 1]
  -v, --version         show program's version number and exit
```

#### Example

```bash
# donwload example dataset
wget https://github.com/akikuno/DAJIN2/raw/main/examples/single.tar.gz
tar -xf single.tar.gz

DAJIN2 \
    --name stx2-deletion \
    --sample "single/barcode25.fq.gz" \
    --control "single/barcode30.fq.gz" \
    --allele "single/design_stx2.fa" \
    --genome mm10 \
    --threads 10

# 2023-06-04 11:30:03: single/barcode30.fq.gz is now processing...
# 2023-06-04 11:30:06: Preprocess single/barcode30.fq.gz...
# 2023-06-04 11:30:06: Mapping single/barcode30.fq.gz...
# 2023-06-04 11:30:21: Call MIDSV single/barcode30.fq.gz...
# 2023-06-04 11:30:31: 🍵 single/barcode30.fq.gz is finished!
# 2023-06-04 11:30:31: single/barcode25.fq.gz is now processing...
# 2023-06-04 11:30:35: Preprocess single/barcode25.fq.gz...
# 2023-06-04 11:34:13: Classify single/barcode25.fq.gz...
# 2023-06-04 11:34:18: Clustering single/barcode25.fq.gz...
# 2023-06-04 11:35:01: Consensus calling single/barcode25.fq.gz...
# 2023-06-04 11:35:08: 🍵 single/barcode25.fq.gz is finished!
# 🎉 Finished! Open DAJINResults/stx2-deletion to see the report.
```

### Batch-mode

DAJIN2 can handle multiple FASTQ files via `batch` subcommand.

```bash
DAJIN2 batch [-h] -f FILE [-t THREADS]

options:
  -h, --help            Show this help message and exit
  -f FILE, --file FILE  CSV or Excel file
  -t THREADS, --threads THREADS
                        Number of threads [default: 1]
```

#### Example

🚧 Working in progress 🚧

## Reference

Kuno A, Ikeda Y, Ayabe S, Kato K, Sakamoto K, et al. (2022) DAJIN enables multiplex genotyping to simultaneously validate intended and unintended target genome editing outcomes. PLOS Biology 20(1): e3001507. https://doi.org/10.1371/journal.pbio.3001507
