[![License](https://img.shields.io/badge/License-MIT-9cf.svg)](https://choosealicense.com/licenses/mit/)
[![Test](https://img.shields.io/github/actions/workflow/status/akikuno/dajin2/pytest.yml?branch=main&label=Test&color=brightgreen)](https://github.com/akikuno/dajin2/actions)
[![Python](https://img.shields.io/pypi/pyversions/DAJIN2.svg?label=Python&color=blue)](https://pypi.org/project/DAJIN2/)
[![PyPI](https://img.shields.io/pypi/v/DAJIN2.svg?label=PyPI&color=orange)](https://pypi.org/project/DAJIN2/)
[![Bioconda](https://img.shields.io/conda/v/bioconda/dajin2?label=Bioconda&color=orange)](https://anaconda.org/bioconda/dajin2)
[![DOI](https://zenodo.org/badge/387721337.svg)](https://zenodo.org/badge/latestdoi/387721337)
[![Paper](https://img.shields.io/badge/Plos%20Biol-10.1371/journal.pbio.3001507-lightgreen)](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3001507)


<p align="center">
<img src="https://user-images.githubusercontent.com/15861316/261833016-7f356960-88cf-4574-87e2-36162b174340.png" width="90%">
</p>

[日本語はこちら](https://github.com/akikuno/DAJIN2/blob/main/docs/README_JP.md)

DAJIN2 is a genotyping tool for genome-edited samples, utilizing nanopore sequencer target sequencing.

The name DAJIN is derived from the phrase 一網**打尽** (Ichimou **DAJIN** or Yīwǎng **Dǎjìn**), symbolizing the concept of capturing everything in one sweep.  

# 🌟 Features

+ **Comprehensive Mutation Detection**: Equipped with the capability to detect genome editing events over a wide range, it can identify a broad spectrum of mutations, from small changes to large structural variations.
  + DAJIN2 is also possible to detect complex mutations characteristic of genome editing, such as "insertions occurring in regions where deletions have occurred."
+ **Intuitive Visualization**: The outcomes of genome editing are visualized intuitively, allowing for the rapid and easy identification and analysis of mutations.
+ **Multi-Sample Compatibility**: Enabling parallel processing of multiple samples. This facilitates efficient progression of large-scale experiments and comparative studies.


# 🛠 Installation

## Prerequisites

- Python >= 3.9
- Unix-like environment (Linux, macOS, WSL2, etc.)

## From [Bioconda](https://anaconda.org/bioconda/DAJIN2) (Recommended)

```bash
# Setup of Bioconda
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority flexible

# Install DAJIN2
conda create -n env-dajin2 python=3.12 DAJIN2 -y
conda activate env-dajin2
```

## From [PyPI](https://pypi.org/project/DAJIN2/)

```bash
pip install DAJIN2
```

> [!IMPORTANT]
> DAJIN2 is actively developed and continuously improved. To access the latest features, please ensure you have the newest version installed.
> ```bash
> DAJIN2 --version
> ```
> To update to the latest version, use one of the following commands:
> ```bash
> conda update DAJIN2 -y
> ```
> or
> ```bash
> pip install -U DAJIN2
> ```


> [!CAUTION]
> If you encounter any issues during the installation, please refer to the [Troubleshooting Guide](https://github.com/akikuno/DAJIN2/blob/main/docs/TROUBLESHOOTING.md)


# 💻 Usage

## Required Files

### FASTQ/FASTA/BAM Files for Sample and Control

In DAJIN2, a **control that has not undergone genome editing** is necessary to detect genome-editing-specific mutations. Specify a directory containing the FASTQ/FASTA (both gzip compressed and uncompressed) or BAM files of the genome editing sample and control.

#### Basecalling with [Guppy](https://community.nanoporetech.com/docs/prepare/library_prep_protocols/Guppy-protocol/v/gpb_2003_v1_revax_14dec2018/guppy-software-overview)
After basecalling with Guppy, the following file structure will be output:

```text
fastq_pass
├── barcode01
│   ├── fastq_runid_b347657c88dced2d15bf90ee6a1112a3ae91c1af_0_0.fastq.gz
│   ├── fastq_runid_b347657c88dced2d15bf90ee6a1112a3ae91c1af_10_0.fastq.gz
│   └── fastq_runid_b347657c88dced2d15bf90ee6a1112a3ae91c1af_11_0.fastq.gz
└── barcode02
    ├── fastq_runid_b347657c88dced2d15bf90ee6a1112a3ae91c1af_0_0.fastq.gz
    ├── fastq_runid_b347657c88dced2d15bf90ee6a1112a3ae91c1af_10_0.fastq.gz
    └── fastq_runid_b347657c88dced2d15bf90ee6a1112a3ae91c1af_11_0.fastq.gz
```

Assuming barcode01 is the control and barcode02 is the sample, the respective directories are specified as follows:

+ Control: `fastq_pass/barcode01`
+ Sample: `fastq_pass/barcode02`


#### Basecalling with [Dorado](https://github.com/nanoporetech/dorado)

For basecalling with Dorado ([`dorado demux`](https://github.com/nanoporetech/dorado?tab=readme-ov-file#barcode-classification)), the following file structure will be output:

```text
dorado_demultiplex
├── EXP-PBC096_barcode01.bam
└── EXP-PBC096_barcode02.bam
```

> [!IMPORTANT]
> Store each BAM file in a separate directory. The directory names can be set arbitrarily.

```text
dorado_demultiplex
├── barcode01
│   └── EXP-PBC096_barcode01.bam
└── barcode02
    └── EXP-PBC096_barcode02.bam
```

Similarly, store the FASTA files outputted after sequence error correction with [`dorado correct`](https://github.com/nanoporetech/dorado) in separate directories.

```text
dorado_correct
├── barcode01
│   └── EXP-PBC096_barcode01.fasta
└── barcode02
    └── EXP-PBC096_barcode02.fasta
```

Assuming barcode01 is the control and barcode02 is the sample, the respective directories are specified as follows:

+ Control: `dorado_demultiplex/barcode01` / `dorado_correct/barcode01`
+ Sample: `dorado_demultiplex/barcode02` / `dorado_correct/barcode02`

### FASTA File Including Anticipated Allele Sequences

The FASTA file should contain descriptions of the alleles anticipated as a result of genome editing.

> [!IMPORTANT]
> **A header name `>control` and its sequence are nessesary.**

If there are anticipated alleles (e.g., knock-ins or knock-outs), include their sequences in the FASTA file too. These anticipated alleles can be named arbitrarily.

Below is an example of a FASTA file:

```text
>control
ACGTACGTACGTACGT
>knock-in
ACGTACGTCCCCACGTACGT
>knock-out
ACGTACGT
```

Here, `>control` represents the sequence of the control allele, while `>knock-in` and `>knock-out` represent the sequences of the anticipated knock-in and knock-out alleles, respectively.

> [!IMPORTANT]
> **Ensure that both ends of the FASTA sequence match those of the amplicon sequence.** If the FASTA sequence is longer or shorter than the amplicon, the difference may be recognized as an indel.  

## Single Sample Analysis

DAJIN2 allows for the analysis of single samples (one sample vs one control).

```bash
DAJIN2 <-s|--sample> <-c|--control> <-a|--allele> <-n|--name> \
  [-g|--genome] [-t|--threads] [-h|--help] [-v|--version]

Options:
-s, --sample              Specify the path to the directory containing sample FASTQ/FASTA/BAM files.
-c, --control             Specify the path to the directory containing control FASTQ/FASTA/BAM files.
-a, --allele              Specify the path to the FASTA file.
-n, --name (Optional)     Set the output directory name. Default: 'Results'.
-g, --genome (Optional)   Specify the reference UCSC genome ID (e.g., hg38, mm39). Default: '' (empty string).
-t, --threads (Optional)  Set the number of threads. Default: 1.
-h, --help                Display this help message and exit.
-v, --version             Display the version number and exit.
```

### Example

```bash
# Download example dataset
curl -LJO https://github.com/akikuno/DAJIN2/raw/main/examples/example_single.tar.gz
tar -xf example_single.tar.gz

# Run DAJIN2
DAJIN2 \
    --control example_single/control \
    --sample example_single/sample \
    --allele example_single/stx2_deletion.fa \
    --name stx2_deletion \
    --genome mm39 \
    --threads 4
```

## Batch Processing

By using the `batch` subcommand, you can process multiple files simultaneously.  
For this purpose, a CSV or Excel file consolidating the sample information is required.  

> [!NOTE]
> For guidance on how to compile sample information, please refer to [this document](https://docs.google.com/presentation/d/e/2PACX-1vSMEmXJPG2TNjfT66XZJRzqJd82aAqO5gJrdEzyhn15YBBr_Li-j5puOgVChYf3jA/embed?start=false&loop=false&delayms=3000).


```bash
DAJIN2 batch <-f|--file> [-t|--threads] [-h]

options:
  -f, --file                Specify the path to the CSV or Excel file.
  -t, --threads (Optional)  Set the number of threads. Default: 1.
  -h, --help                Display this help message and exit.
```

### Example

```bash
# Donwload the example dataset
curl -LJO https://github.com/akikuno/DAJIN2/raw/main/examples/example_batch.tar.gz
tar -xf example_batch.tar.gz

# Run DAJIN2
DAJIN2 batch --file example_batch/batch.csv --threads 4
```


# 📈 Reports

Upon completion of DAJIN2 processing, a directory named `DAJIN_Results/{NAME}` is generated.  
Inside the `DAJIN_Results/{NAME}` directory, the following files can be found:  

```
DAJIN_Results/tyr-substitution
├── BAM
│   ├── tyr_c230gt_01
│   ├── tyr_c230gt_10
│   ├── tyr_c230gt_50
│   └── tyr_control
├── FASTA
│   ├── tyr_c230gt_01
│   ├── tyr_c230gt_10
│   └── tyr_c230gt_50
├── HTML
│   ├── tyr_c230gt_01
│   ├── tyr_c230gt_10
│   └── tyr_c230gt_50
├── MUTATION_INFO
│   ├── tyr_c230gt_01.csv
│   ├── tyr_c230gt_10.csv
│   └── tyr_c230gt_50.csv
├── read_plot.html
├── read_plot.pdf
└── read_summary.xlsx
```

## 1. BAM

The BAM directory contains the BAM files of reads classified per allele.  

> [!NOTE]
> Specifying a reference genome using the `genome` option will align the reads to that genome.  
> Without `genome` options, the reads will align to the control allele within the input FASTA file.

## 2. FASTA and HTML

The FASTA directory stores the FASTA files of each allele.  
The HTML directory contains HTML files for each allele, where mutation sites are color-highlighted.  
For example, Tyr point mutation is highlighted in **green**.  

<img src="https://raw.githubusercontent.com/akikuno/logos/refs/heads/main/DAJIN2/tyr-substitution.png" width="75%" />

Furthermore, DAJIN2 extracts representative SV alleles (Insertion, Deletion, Inversion) included in the sample and highlights SV regions with colored underlines.  
The following is an example where a deletion (light blue) and an insertion (red) are observed at both ends of an inversion (purple underline):

<img src="https://raw.githubusercontent.com/akikuno/logos/refs/heads/main/DAJIN2/cables2-inversion.png" width="75%" />

## 3. MUTATION_INFO

The MUTATION_INFO directory saves tables depicting mutation sites for each allele.  
An example of a *Tyr* point mutation is described by its position on the chromosome and the type of mutation.  

<img src="https://user-images.githubusercontent.com/15861316/274519342-a613490d-5dbb-4a27-a2cf-bca0686b30f0.png" width="75%">

## 4. resd_summary.xlsx, read_plot.html and read_plot.pdf

read_summary.xlsx describes the number of reads and presence proportion for each allele.  
Both read_plot.html and read_plot.pdf illustrate the proportions of each allele.  
The chart's **Allele type** indicates the type of allele, and **Percent of reads** shows the proportion of reads for each allele.  

The **Allele type** includes:
- **Intact**: Alleles that perfectly match the input FASTA allele.
- **Indels**: Substitutions, deletions, insertions, or inversions within 50 bases.
- **SV**: Substitutions, deletions, insertions, or inversions beyond 50 bases.

<img src="https://user-images.githubusercontent.com/15861316/274521067-4d217251-4c62-4dc9-9c05-7f5377dd3025.png" width="75%">

> [!WARNING]
> In PCR amplicon sequencing, the % of reads might not match the actual allele proportions due to amplification bias.  
> Especially when large deletions are present, the deletion alleles might be significantly amplified, potentially not reflecting the actual allele proportions.

# 📣 Feedback and Support

> [!NOTE]
> For frequently asked questions, please refer to [this page](https://github.com/akikuno/DAJIN2/blob/main/docs/FAQ.md).


For more questions, bug reports, or other forms of feedback, we'd love to hear from you!  
Please use [GitHub Issues](https://github.com/akikuno/DAJIN2/issues/new/choose) for all reporting purposes.  

Please refer to [CONTRIBUTING](https://github.com/akikuno/DAJIN2/blob/main/docs/CONTRIBUTING.md) for how to contribute and how to verify your contributions.  


# 🤝 Code of Conduct

Please note that this project is released with a [Contributor Code of Conduct](https://github.com/akikuno/DAJIN2/blob/main/docs/CODE_OF_CONDUCT.md).  
By participating in this project you agree to abide by its terms.  


# 📄 References

For more information, please refer to the following publication:

[Kuno A, et al. (2022) DAJIN enables multiplex genotyping to simultaneously validate intended and unintended target genome editing outcomes. *PLoS Biology* 20(1): e3001507.](https://doi.org/10.1371/journal.pbio.3001507)

