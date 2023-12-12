[![License](https://img.shields.io/badge/License-MIT-9cf.svg?style=flat-square)](https://choosealicense.com/licenses/mit/)
[![Test](https://img.shields.io/github/actions/workflow/status/akikuno/dajin2/pytest.yml?branch=main&label=Test&color=brightgreen&style=flat-square)](https://github.com/akikuno/dajin2/actions)
[![Python](https://img.shields.io/pypi/pyversions/DAJIN2.svg?label=Python&color=blue&style=flat-square)](https://pypi.org/project/DAJIN2/)
[![PyPI](https://img.shields.io/pypi/v/DAJIN2.svg?label=PyPI&color=orange&style=flat-square)](https://pypi.org/project/DAJIN2/)
[![Bioconda](https://img.shields.io/conda/v/bioconda/dajin2?label=Bioconda&color=orange&style=flat-square)](https://anaconda.org/bioconda/dajin2)
[![DOI](https://zenodo.org/badge/387721337.svg)](https://zenodo.org/badge/latestdoi/387721337)


<p align="center">
<img src="https://user-images.githubusercontent.com/15861316/261833016-7f356960-88cf-4574-87e2-36162b174340.png" width="90%">
</p>

[日本語はこちら](https://github.com/akikuno/DAJIN2/blob/main/docs/README_JP.md)

DAJIN2 is a genotyping software designed for organisms that have undergone genome editing, utilizing nanopore sequencing technology.  

The name DAJIN is inspired by the term 一網**打尽** (Ichimou **DAJIN** or Yīwǎng **Dǎjìn**), which signifies capturing everything in a single net.  


## 🛠 Installation

### From [Bioconda](https://anaconda.org/bioconda/DAJIN2) (Recommended)

```bash
conda install -c bioconda DAJIN2
```

### From [PyPI](https://pypi.org/project/DAJIN2/)

```bash
pip install DAJIN2
```

> [!CAUTION]
> If you encounter any issues during the installation, please refer to the [Troubleshooting Guide](https://github.com/akikuno/DAJIN2/blob/main/docs/TROUBLESHOOTING.md)


## 💡 Usage

### Single Sample Analysis

DAJIN2 allows for the analysis of single samples (one sample vs one control).

```bash
DAJIN2 <-s|--sample> <-c|--control> <-a|--allele> <-n|--name> \
  [-g|--genome] [-t|--threads] [-h|--help] [-v|--version]

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
# 🎉 Finished! Open DAJIN_Results/stx2-deletion to see the report.
```

### Batch Processing

By using the `batch` subcommand, you can process multiple FASTQ files simultaneously.  
For this purpose, a CSV or Excel file consolidating the sample information is required.  
For a specific example, please refer to [this link](https://github.com/akikuno/DAJIN2/blob/main/examples/example-batch/batch.csv).


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
# 🎉 Finished! Open DAJIN_Results/tyr-substitution to see the report.
```

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


## 📄 References

For more information, please refer to the following publication:

[Kuno A, et al. (2022) DAJIN enables multiplex genotyping to simultaneously validate intended and unintended target genome editing outcomes. *PLoS Biology* 20(1): e3001507.](https://doi.org/10.1371/journal.pbio.3001507)


## 📣Feedback and Support

For questions, bug reports, or other forms of feedback, we'd love to hear from you!  
Please use [GitHub Issues](https://github.com/akikuno/DAJIN2/issues) for all reporting purposes.  

Please refer to [CONTRIBUTING](https://github.com/akikuno/DAJIN2/blob/main/docs/CONTRIBUTING.md) for how to contribute and how to verify your contributions.  

## 🤝 Code of Conduct

Please note that this project is released with a [Contributor Code of Conduct](https://github.com/akikuno/DAJIN2/blob/main/docs/CODE_OF_CONDUCT.md).  
By participating in this project you agree to abide by its terms.  
