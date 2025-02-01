# Troubleshooting Guide

## Installation

### Long Story Short

Please run the following command. If an error occurs, please report it on [GitHub Issues](https://github.com/akikuno/DAJIN2/issues/new?template=bug_report.yml).

```bash
# Update conda
conda update -n base conda -y

# Setup of Bioconda
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority flexible # <- IMPORTANT!!

# Install DAJIN2 to a virtual environment
conda create -n env-dajin2 python=3.12 -y
conda activate env-dajin2
conda install -c conda-forge -c bioconda DAJIN2 -y
```

以下の項目は、詳細です。

### Prerequisites

Before installing `DAJIN2`, please ensure that your system meets the following requirements:

- Maintained Python versions (excluding the latest version)
- Unix-like environment (Linux, macOS, WSL2, etc.)
- [Conda](https://docs.conda.io/en/latest/) is highly recommended for managing dependencies
- If using pip, access to administrative privileges or the ability to install packages globally

### Dependencies

`DAJIN2` depends on the `pysam` package, which in turn requires `htslib`. These dependencies are critical for the functionality of `DAJIN2` but can pose installation challenges:

- `htslib` requires `zlib.h`, which may not be available on all systems by default.
- As of the latest update, `htslib v1.21`, [there is no support for Windows via Bioconda](https://anaconda.org/bioconda/htslib).

### Recommended Installation Strategy

#### Using Conda

We strongly recommend using Conda for installation, as they efficiently manage the complex dependencies of `pysam` and `htslib`:

1. Install the latest Conda if you have not already. You can download `miniforge` from [here](https://github.com/conda-forge/miniforge?tab=readme-ov-file#miniforge3).

2. Setup the [Bioconda](https://bioconda.github.io/) channel:

```bash
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority flexible
```

> [!CAUTION]
> Bioconda recommends using `conda config --set channel_priority strict`, but as the installation of openpyxl fails, please use **`conda config --set channel_priority flexible`**.

> [!WARNING]
> If you encounter the following error message when you install DAJIN2, please try **`conda config --set channel_priority flexible`**.
> ```bash
> LibMambaUnsatisfiableError: Encountered problems while solving:
> package dajin2-x.x.x-xxx is excluded by strict repo priority
> ```

3. Create a new virtual environment (optional but recommended):

```bash
conda create -n env-dajin2
conda activate env-dajin2
```

4. Install `DAJIN2` in the virtual environment:
```bash
conda install -c bioconda DAJIN2
```

#### Using pip

If you prefer or are required to use pip, please ensure that `zlib.h` and other dependencies are properly installed on your system. This method might require administrative privileges:

```bash
pip install DAJIN2
```

> [!CAUTION]
> Pip installation might encounter issues due to the dependencies mentioned above, especially on systems without `zlib.h` or on Windows.
> ```bash
> sudo apt install gcc zlib1g zlib1g-dev # WSL2 or Ubuntu  
> brew install gcc zlib # macOS
> ```


## Report other troubles

Please use [GitHub Issues](https://github.com/akikuno/DAJIN2/issues) for all reporting purposes.  
