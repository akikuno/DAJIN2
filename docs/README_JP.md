# DAJIN2

⚠️ DAJIN2は現在開発中です。2023年3月のリリースを目指しています。

## 使い方

```bash
DAJIN2 [-s FASTQ] [-c FASTQ] [-a FASTA] [-o OUTPUT] [-g GENOME] [-t THREADS]
```

### Batchモード

```bash
DAJIN2 batch [-f <batch.csv>] [-t THREADS]
```

### 引数

```
-s, --sample    : (required) Give the full path to FASTQ file(s) (uncompressed or gzip compressed)

-c, --control   : (required) Give the full path to FASTQ file(s) (uncompressed or gzip compressed)

-a, --allele    : (required) Give the full path to FASTA file(s). `>control` sequence is required. (uncompressed)

-o, --output    : (required) Give the full path to the output vcf file [default: DAJIN_results ]

-g, --genome    : (optional) Specify the reference genome listed in UCSC Genome Browser Gateway (http://hgdownload.soe.ucsc.edu/downloads.html)

-t, --threads   : (optional) Maximum number of threads to use [default: 1]

-h, --help      : (optional) Show the help message and exit

-v, --version   : (optional) Show the version and exit
```
