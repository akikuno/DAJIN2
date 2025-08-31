[![License](https://img.shields.io/badge/License-MIT-9cf.svg)](https://choosealicense.com/licenses/mit/)
[![Test](https://img.shields.io/github/actions/workflow/status/akikuno/dajin2/pytest.yml?branch=main&label=Test&color=brightgreen)](https://github.com/akikuno/dajin2/actions)
[![Python](https://img.shields.io/pypi/pyversions/DAJIN2.svg?label=Python&color=blue)](https://pypi.org/project/DAJIN2/)
[![PyPI](https://img.shields.io/pypi/v/DAJIN2.svg?label=PyPI&color=orange)](https://pypi.org/project/DAJIN2/)
[![Bioconda](https://img.shields.io/conda/v/bioconda/dajin2?label=Bioconda&color=orange)](https://anaconda.org/bioconda/dajin2)
[![DOI](https://zenodo.org/badge/387721337.svg)](https://zenodo.org/badge/latestdoi/387721337)
[![Paper](https://img.shields.io/badge/Plos%20Biol-10.1371/journal.pbio.3001507-lightgreen)](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3001507)
[![Contact💬](https://img.shields.io/badge/Contact💬-923DE2)](https://forms.gle/r4YRs1th7NGHfDcS9)


<p align="center">
<img src="https://user-images.githubusercontent.com/15861316/261833016-7f356960-88cf-4574-87e2-36162b174340.png" width="90%">
</p>

[日本語版 README はこちら](https://github.com/akikuno/DAJIN2/blob/main/docs/README_JP.md)

DAJIN2 is a genotyping tool for genome-edited samples, utilizing nanopore target sequencing.

**DAJIN2** takes its name from the Japanese phrase 一網**打尽** (*Ichimou DAJIN*, or *Yīwǎng Dǎjìn* in Chinese),  
which means “to capture everything in a single sweep.”  
This reflects the tool’s design philosophy: to comprehensively detect both intended and unintended genome editing outcomes in one go.

# 🌟 Features

+ **Comprehensive Mutation Detection**  
  DAJIN2 can detect a wide range of genome editing events in nanopore-targeted regions, from point mutations to structural variants.  
  It is particularly effective at identifying **unexpected mutations** and **complex mutations**, such as insertions within deleted regions.

+ **Highly Sensitive Allele Classification**  
  Supports classification of mosaic alleles, capable of detecting minor alleles present at approximately 1%.

+ **Intuitive Visualization**  
  Genome editing results are visualized in an intuitive manner, enabling rapid and easy identification of mutations.

+ **Multi-Sample Support**  
  Batch processing of multiple samples is supported, allowing efficient execution of large-scale experiments and comparative studies.

+ **Simple Installation and Operation**  
  Requires no specialized computing environment and runs smoothly on a standard laptop.  
  Easily installable via Bioconda or PyPI, and usable via the command line.  

# 🛠 Installation

## System Requirements

### Hardware

- **Runs on a standard laptop**
- Recommended memory: 8 GB or more

>[!NOTE]
> DAJIN2 is the successor to DAJIN, which required a GPU for efficient computation due to its use of deep learning.  
> In contrast, **DAJIN2 does not use deep learning and does not require a GPU**.  
> Therefore, it runs smoothly on typical laptops.

### Software

- Python 3.9-3.12
- Unix-based environment (Linux, macOS, WSL2, etc.)

>[!IMPORTANT]
> **For Windows Users**  
> DAJIN2 is designed to run in a Linux environment.  
> If you are using Windows, please use **WSL2 (Windows Subsystem for Linux 2)**.


## From [Bioconda](https://anaconda.org/bioconda/DAJIN2) (Recommended)

```bash
# Setting up Bioconda
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
> **DAJIN2 is actively being developed and improved.**  
> Please make sure you are using the latest version to take advantage of the newest features.  
>
> 🔍 **To check your current version:**
> ```bash
> DAJIN2 --version
> ```
>
> ➡️ **Check the latest version:**  
> https://github.com/akikuno/DAJIN2/releases
>
> 🔄 **To update to the latest version:**
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
> **A header name `>control` and its sequence are necessary.**

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
  [-g|--genome] [-b|--bed] [-t|--threads] [--no-filter] [-h|--help] [-v|--version]

Options:
-s, --sample            Specify the path to the directory containing sample FASTQ/FASTA/BAM files.
-c, --control           Specify the path to the directory containing control FASTQ/FASTA/BAM files.
-a, --allele            Specify the path to the FASTA file.
-n, --name (Optional)   Set the output directory name. Default: 'Results'.
-g, --genome (Optional) Specify the reference UCSC genome ID (e.g., hg38, mm39). Default: '' (empty string).
-b, --bed (Optional)    Specify the path to BED6 file containing genomic coordinates. Default: '' (empty string).
-t, --threads (Optional) Set the number of threads. Default: 1.
--no-filter (Optional)  Disable minor allele filtering (keep alleles <0.5%). Default: False.
-h, --help              Display this help message and exit.
-v, --version           Display the version number and exit.
```

### Using BED Files for Genomic Coordinates

If the reference genome is not from UCSC, or if the external servers that DAJIN2 depends on (UCSC Genome Browser and GGGENOME) are unavailable, you can specify a BED file using the `-b/--bed` option to run offline.

When using the `-b/--bed` option with a BED file, please ensure:

1. **Use BED6 format** (6 columns required):
   ```
   chr1    1000000    1001000    mm39    248956422    +
   ```
   
   **Column descriptions:**
   - Column 1: Chromosome name (e.g., chr1, chr2)
   - Column 2: Start position (0-based, inclusive)
   - Column 3: End position (0-based, exclusive)
   - Column 4: Name (**genome ID**)
   - Column 5: Score (**chromosome size for proper IGV visualization**)
   - Column 6: Strand (+ or -, **must match FASTA allele orientation**)

> [!NOTE]  
> For the score field (column 5), please enter the size of the chromosome specified in column 1.  
> While the original BED format limits scores to 1000, DAJIN2 accepts **chromosome sizes without any issue**.

> [!IMPORTANT]  
> **Strand orientation must match**. The strand field (column 6: `+` or `-`) in your BED file **must match the strand orientation of your FASTA allele sequences**.  
> - If your FASTA allele sequence is on the **forward strand** (5' to 3'), use `+` in the BED file  
> - If your FASTA allele sequence is on the **reverse strand** (3' to 5'), use `-` in the BED file

For detailed BED file usage, see [BED_COORDINATE_USAGE.md](docs/BED_COORDINATE_USAGE.md).

### Rare Mutation Detection with `--no-filter`

By default, DAJIN2 filters out alleles with read counts below 0.5% (5 reads out of 100,000 downsampled reads) to reduce noise and improve accuracy. However, when analyzing rare mutations or somatic mosaicism where minor alleles may be present at very low frequencies, you can use the `--no-filter` option to disable this filtering.

**When to use `--no-filter`:**
- Detecting rare somatic mutations (< 0.5% frequency)
- Analyzing samples with suspected low-level mosaicism
- Research requiring detection of all possible alleles regardless of frequency

**Usage:**
```bash
DAJIN2 \
    --control example_single/control \
    --sample example_single/sample \
    --allele example_single/stx2_deletion.fa \
    --name stx2_deletion \
    --genome mm39 \
    --threads 4 \
    --no-filter
```

> [!CAUTION]
> Using `--no-filter` may increase noise and false positives in the results. It is recommended to validate rare alleles through additional experimental methods.

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

By using the `batch` subcommand, you can process multiple samples simultaneously.  
For this purpose, a CSV or Excel file consolidating the sample information is required.  

> [!NOTE]
> For guidance on how to compile sample information, please refer to [this document](https://docs.google.com/presentation/d/e/2PACX-1vSMEmXJPG2TNjfT66XZJRzqJd82aAqO5gJrdEzyhn15YBBr_Li-j5puOgVChYf3jA/embed?start=false&loop=false&delayms=3000).

**Required columns:** `sample`, `control`, `allele`, `name`  
**Optional columns:** `genome`, `bed` (or `genome_coordinate`), and any custom columns

**Example CSV with BED files:**
```csv
sample,control,allele,name,bed
/path/to/sample1,/path/to/control1,/path/to/allele1.fa,experiment1,/path/to/coords1.bed
/path/to/sample2,/path/to/control2,/path/to/allele2.fa,experiment2,/path/to/coords2.bed
```

> [!TIP]
> It is **recommended to use the same value in the `name` column for samples that belong to the same experiment.**  
> Using identical names enables parallel processing, thereby improving efficiency.  
> Here's an example 👉 [batch.csv](https://github.com/akikuno/DAJIN2/blob/main/examples/example_batch/batch.csv)


```bash
DAJIN2 batch <-f|--file> [-t|--threads] [--no-filter] [-h]

options:
  -f, --file                Specify the path to the CSV or Excel file.
  -t, --threads (Optional)  Set the number of threads. Default: 1.
  --no-filter (Optional)    Disable minor allele filtering (keep alleles <0.5%). Default: False.
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

## GUI (Graphical User Interface) Mode

DAJIN2 provides a web interface that can be launched with a single command:

```bash
DAJIN2 gui
```

When executed, your default web browser will open and display the following GUI at `http://localhost:{PORT}/`.

<img src="https://raw.githubusercontent.com/akikuno/DAJIN2/refs/heads/main/image/dajin2-gui.jpg" width="75%">

> [!NOTE]
> If the browser does not launch automatically, please open your browser manually and navigate to `http://localhost:{PORT}/`.

### Single Sample Analysis via GUI

1. **Launch GUI**  
   Run `DAJIN2 gui` to open the web interface.

2. **Project Setup**  
   - **Project Name**: Enter any analysis name  
   - **Directory Upload**: Select directories containing sample or control FASTQ/FASTA/BAM files  
   - **Allele FASTA**: Upload FASTA file containing expected allele sequences  
   - **BED File (optional)**: Upload BED6 format file to specify genomic coordinates

3. **Parameter Configuration**  
   - **Reference Genome (optional)**: Specify UCSC genome ID (e.g., `hg38`, `mm39`)  
   - **Threads**: Set the number of CPU threads to use  
   - **No Filter**: Enable to detect rare mutations below 0.5% frequency

4. **Run Analysis**  
   Click "Start Analysis" and the progress will be displayed in real-time.

5. **View Results**  
   After completion, the output folder path will be displayed for accessing result files.

### Batch Processing via GUI

1. **Prepare Batch File**  
   Create a CSV or Excel file with columns: `sample`, `control`, `allele`, `name`.

2. **Upload Batch File**  
   Use the "Batch Processing" tab to upload your configuration file.

3. **Configure Global Settings**  
   Set threads and filtering options for all samples at once.

4. **Monitor Progress**  
   The analysis status for each sample is displayed with detailed log output.

5. **View Results**  
   Results are saved in the `DAJIN_Results/` folder with subdirectories for each sample.


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

The BAM directory contains the BAM files of reads classified by allele.  

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

## 4. read_summary.xlsx, read_plot.html and read_plot.pdf

read_summary.xlsx summarizes the number and proportion of reads per allele.  
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

We welcome your questions, bug reports, and feedback.  
Please use the following Google Form to submit your report:  
👉 [Google Form](https://forms.gle/r4YRs1th7NGHfDcS9)  

If you have a GitHub account, you can also submit reports via  
👉 [GitHub Issues](https://github.com/akikuno/DAJIN2/issues/new/choose)  

Please refer to [CONTRIBUTING](https://github.com/akikuno/DAJIN2/blob/main/docs/CONTRIBUTING.md) for how to contribute and how to verify your contributions.  

> [!NOTE]
> For frequently asked questions, please refer to [this page](https://github.com/akikuno/DAJIN2/blob/main/docs/FAQ.md).


# 🤝 Code of Conduct

Please note that this project is released with a [Contributor Code of Conduct](https://github.com/akikuno/DAJIN2/blob/main/docs/CODE_OF_CONDUCT.md).  
By participating in this project you agree to abide by its terms.  


# 📄 References

For more information, please refer to the following publication:

[Kuno A, et al. (2022) DAJIN enables multiplex genotyping to simultaneously validate intended and unintended target genome editing outcomes. *PLoS Biology* 20(1): e3001507.](https://doi.org/10.1371/journal.pbio.3001507)

