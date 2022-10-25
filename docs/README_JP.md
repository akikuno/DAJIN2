# DAJIN2

⚠️ DAJIN2は現在開発中です。2023年3月のリリースを予定しています。⚠️

# 使い方

## Singleモード

```bash
DAJIN2 \
    [-s/--sample] FASTQ \
    [-c/--control] FASTQ \
    [-a/--allele] FASTA \
    [-n/--name] NAME \
    [-g/--genome] GENOME \
    [-t/--threads] THREADS
```


### 引数

```
-s, --sample    : (required) Give the full path to FASTQ file(s) (uncompressed or gzip compressed)

-c, --control   : (required) Give the full path to a FASTQ file (uncompressed or gzip compressed)

-a, --allele    : (required) Give the full path to a FASTA file. `>control` sequence is required. (uncompressed)

-o, --output    : (required) Give the full path to the output vcf file [default: DAJIN_results ]

-g, --genome    : (optional) Specify the reference genome listed in UCSC Genome Browser Gateway (http://hgdownload.soe.ucsc.edu/downloads.html)

-t, --threads   : (optional) Maximum number of threads to use [default: 1]

-h, --help      : (optional) Show the help message and exit

-v, --version   : (optional) Show the version and exit
```

## Batchモード

```bash
DAJIN2 batch [-f/--file] batch.csv [-t/--threads] THREADS
```
