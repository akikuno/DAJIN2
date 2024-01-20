[![License](https://img.shields.io/badge/License-MIT-9cf.svg)](https://choosealicense.com/licenses/mit/)
[![Test](https://img.shields.io/github/actions/workflow/status/akikuno/dajin2/pytest.yml?branch=main&label=Test&color=brightgreen)](https://github.com/akikuno/dajin2/actions)
[![Python](https://img.shields.io/pypi/pyversions/DAJIN2.svg?label=Python&color=blue)](https://pypi.org/project/DAJIN2/)
[![PyPI](https://img.shields.io/pypi/v/DAJIN2.svg?label=PyPI&color=orange)](https://pypi.org/project/DAJIN2/)
[![Bioconda](https://img.shields.io/conda/v/bioconda/dajin2?label=Bioconda&color=orange)](https://anaconda.org/bioconda/dajin2)
[![DOI](https://zenodo.org/badge/387721337.svg)](https://zenodo.org/badge/latestdoi/387721337)


<p align="center">
<img src="https://user-images.githubusercontent.com/15861316/261833016-7f356960-88cf-4574-87e2-36162b174340.png" width="90%">
</p>

[日本語はこちら](https://github.com/akikuno/DAJIN2/blob/main/docs/README_JP.md)

DAJIN2 is a genotyping tool for genome-edited samples, utilizing nanopore sequencer target sequencing.

The name DAJIN is derived from the phrase 一網**打尽** (Ichimou **DAJIN** or Yīwǎng **Dǎjìn**), symbolizing the concept of capturing everything in one sweep.  

## 🌟 Features

+ **Comprehensive Mutation Detection**: Equipped with the capability to detect genome editing events over a wide range, it can identify a broad spectrum of mutations, from small changes to large structural variations.
+ **Intuitive Visualization**: The outcomes of genome editing are visualized intuitively, allowing for the rapid and easy identification and analysis of mutations.
+ **Multi-Sample Compatibility**: Accommodates a variety of samples, enabling simultaneous processing of multiple samples. This facilitates efficient progression of large-scale experiments and comparative studies.


## 🛠 Installation

### Prerequisites

- Python 3.7 or later
- Unix-like environment (Linux, macOS, WSL2, etc.)

### From [Bioconda](https://anaconda.org/bioconda/DAJIN2) (Recommended)

```bash
conda create -n env-dajin2 -c conda-forge -c bioconda python=3.10 DAJIN2 -y
conda activate env-dajin2
```

> [!NOTE]
> To Apple Silicon (ARM64) users:  
> [Since the Bioconda channel does not yet support Apple Silicon](https://github.com/bioconda/bioconda-recipes/issues/37068#issuecomment-1257790919), please use the following command to install `DAJIN2` through Rosetta.
> ```bash
> CONDA_SUBDIR=osx-64 conda create -n env-dajin2 -c conda-forge -c bioconda python=3.10 DAJIN2 -y
> conda activate env-dajin2
> conda config --env --set subdir osx-64
> ```

### From [PyPI](https://pypi.org/project/DAJIN2/)

```bash
pip install DAJIN2
```

> [!CAUTION]
> If you encounter any issues during the installation, please refer to the [Troubleshooting Guide](https://github.com/akikuno/DAJIN2/blob/main/docs/TROUBLESHOOTING.md)


## 💡 Usage

### Required Files

#### FASTQ Files for Sample and Control

In DAJIN2, a **control that has not undergone genome editing** is necessary to detect genome-editing-specific mutations. Specify a directory containing the FASTQ files (both gzip compressed and uncompressed) of the genome editing sample and control.

<!-- [Nanopore Guppy](https://community.nanoporetech.com/docs/prepare/library_prep_protocols/Guppy-protocol) -->
After base calling with Guppy, the following file structure will be output:


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

Assuming barcode01 as the control and barcode02 as the sample, specify each directory as follows:

+ Control: fastq_pass/barcode01
+ Sample: fastq_pass/barcode01

#### FASTA File for Allele

The FASTA file should contain descriptions of the alleles anticipated as a result of genome editing.

> [!IMPORTANT]
> Specifying the control allele: A header name >control and its sequence are mandatory.

If there are anticipated alleles (e.g., knock-ins or knock-outs), include their sequences in the FASTA file too. These anticipated alleles can be named arbitrarily.

Below is a typical example of a FASTA file:

```text
>control
ACGTACGTACGTACGT
>knock-in
ACGTACGTCCCCACGTACGT
>knock-out
ACGTACGT
```

Here, `>control` represents the sequence of the control allele, while `>knock-in` and `>knock-out` represent the sequences of the anticipated knock-in and knock-out alleles, respectively.


### Single Sample Analysis

DAJIN2 allows for the analysis of single samples (one sample vs one control).

```bash
DAJIN2 <-s|--sample> <-c|--control> <-a|--allele> <-n|--name> \
  [-g|--genome] [-t|--threads] [-h|--help] [-v|--version]

Options:
-s, --sample              Specify the path to the directory containing sample FASTQ files.
-c, --control             Specify the path to the directory containing control FASTQ files.
-a, --allele              Specify the path to the FASTA file.
-n, --name (Optional)     Set the output directory name. Default: 'Results'.
-g, --genome (Optional)   Specify the reference genome ID (e.g., hg38, mm39). Default: '' (empty string).
-t, --threads (Optional)  Set the number of threads. Default: 1.
-h, --help                Display this help message and exit.
-v, --version             Display the version number and exit.
```

#### Example

```bash
DAJIN2 \
    --control example/barcode01 \
    --sample example/barcode02 \
    --allele example/design.fa \
    --name IL6-knockin \
    --genome hg38 \
    --threads 4
```

<!-- ```bash
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
# 🎉 Finished! Open DAJIN_Results/stx2-deletion to see the report.
```
 -->

### Batch Processing

By using the `batch` subcommand, you can process multiple FASTQ files simultaneously.  
For this purpose, a CSV or Excel file consolidating the sample information is required.  
<!-- For a specific example, please refer to [this link](https://github.com/akikuno/DAJIN2/blob/main/examples/example-batch/batch.csv). -->

> [!NOTE]
> For guidance on how to compile sample information, please refer to [this document](https://docs.google.com/presentation/d/e/2PACX-1vSMEmXJPG2TNjfT66XZJRzqJd82aAqO5gJrdEzyhn15YBBr_Li-j5puOgVChYf3jA/embed?start=false&loop=false&delayms=3000).


```bash
DAJIN2 batch <-f|--file> [-t|--threads] [-h]

options:
  -f, --file                Specify the path to the CSV or Excel file.
  -t, --threads (Optional)  Set the number of threads. Default: 1.
  -h, --help                Display this help message and exit.
```

#### Example

```bash
DAJIN2 --file batch.csv --threads 4
```

<!-- ```bash
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
# 🎉 Finished! Open DAJIN_Results/tyr-substitution to see the report.
``` -->

## 📈 Report Contents

Upon completion of DAJIN2 processing, a directory named **DAJIN_Results** is generated.  
Inside the **DAJIN_Results** directory, the following files can be found:  

```
DAJIN_Results/tyr-substitution
├── BAM
│   ├── tyr_c230gt_01%
│   ├── tyr_c230gt_10%
│   ├── tyr_c230gt_50%
│   └── tyr_control
├── FASTA
│   ├── tyr_c230gt_01%
│   ├── tyr_c230gt_10%
│   └── tyr_c230gt_50%
├── HTML
│   ├── tyr_c230gt_01%
│   ├── tyr_c230gt_10%
│   └── tyr_c230gt_50%
├── MUTATION_INFO
│   ├── tyr_c230gt_01%.csv
│   ├── tyr_c230gt_10%.csv
│   └── tyr_c230gt_50%.csv
├── read_all.csv
├── read_plot.html
├── read_plot.pdf
└── read_summary.csv
```

### 1. BAM

The BAM directory contains the BAM files of reads classified per allele.  

> [!NOTE]
> Specifying a reference genome using the `genome` option will align the reads to that genome.  
> Without `genome` options, the reads will align to the control allele within the input FASTA file.

### 2. FASTA and HTML

The FASTA directory stores the FASTA files of each allele.  
The HTML directory contains HTML files for each allele, where mutation sites are color-highlighted.  
For example, Tyr point mutation is highlighted in **green**.  

<img src="https://user-images.githubusercontent.com/15861316/274518501-2ca3f442-1b86-4635-be3d-fd37575c4ca2.png" width="75%" />

### 3. MUTATION_INFO

The MUTATION_INFO directory saves tables depicting mutation sites for each allele.  
An example of a Tyr point mutation is described by its position on the chromosome and the type of mutation.  

<img src="https://user-images.githubusercontent.com/15861316/274519342-a613490d-5dbb-4a27-a2cf-bca0686b30f0.png" width="75%">

### 4. read_plot.html and read_plot.pdf

Both read_plot.html and read_plot.pdf illustrate the proportions of each allele.  
The chart's **Allele type** indicates the type of allele, and **% of reads** shows the proportion of reads for that allele.  

Additionally, the types of **Allele type** include:
- **intact**: Alleles that perfectly match the input FASTA allele.
- **indels**: Substitutions, deletions, insertions, or inversions within 50 bases.
- **sv**: Substitutions, deletions, insertions, or inversions beyond 50 bases.

<img src="https://user-images.githubusercontent.com/15861316/274521067-4d217251-4c62-4dc9-9c05-7f5377dd3025.png" width="75%">

> [!WARNING]
> In PCR amplicon sequencing, the % of reads might not match the actual allele proportions due to amplification bias.  
> Especially when large deletions are present, the deletion alleles might be significantly amplified, potentially not reflecting the actual allele proportions.

### 5. read_all.csv and read_summary.csv

- read_all.csv: Records which allele each read is classified under.  
- read_summary.csv: Describes the number of reads and presence proportion for each allele.  

## 📣Feedback and Support

For questions, bug reports, or other forms of feedback, we'd love to hear from you!  
Please use [GitHub Issues](https://github.com/akikuno/DAJIN2/issues) for all reporting purposes.  

Please refer to [CONTRIBUTING](https://github.com/akikuno/DAJIN2/blob/main/docs/CONTRIBUTING.md) for how to contribute and how to verify your contributions.  

## 🤝 Code of Conduct

Please note that this project is released with a [Contributor Code of Conduct](https://github.com/akikuno/DAJIN2/blob/main/docs/CODE_OF_CONDUCT.md).  
By participating in this project you agree to abide by its terms.  


## 📄 References

For more information, please refer to the following publication:

[Kuno A, et al. (2022) DAJIN enables multiplex genotyping to simultaneously validate intended and unintended target genome editing outcomes. *PLoS Biology* 20(1): e3001507.](https://doi.org/10.1371/journal.pbio.3001507)

