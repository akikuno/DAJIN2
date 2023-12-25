# Troubleshooting Guide

## Installation

### Long Story Short

```bash
# Update conda
conda update -n base conda -y

# Setup of Bioconda
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict

# Install DAJIN2 to a virtual environment
conda create -n env-dajin2 -y
conda activate env-dajin2
conda install -c bioconda DAJIN2 -y
```

### Prerequisites

Before installing `DAJIN2`, please ensure that your system meets the following requirements:

- Python >=3.7
- Unix-like environment (Linux, macOS, WSL, etc.)
- [Conda](https://docs.conda.io/en/latest/) is highly recommended for managing dependencies
- If using pip, access to administrative privileges or the ability to install packages globally

### Dependencies

`DAJIN2` depends on the `pysam` package, which in turn requires `htslib`. These dependencies are critical for the functionality of `DAJIN2` but can pose installation challenges:

- `htslib` requires `zlib.h`, which may not be available on all systems by default.
- As of the latest update, `htslib v1.18`, there is no support for Windows via Bioconda.

### Recommended Installation Method

#### Using Conda/Mamba

We strongly recommend using Conda or Mamba for installation, as they efficiently manage the complex dependencies of `pysam` and `htslib`:

1. Install the latest Conda if you have not already. You can download Conda from [here](https://docs.conda.io/en/latest/miniconda.html).

2. Setup the [Bioconda](https://bioconda.github.io/) channel:

```bash
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
```

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

> [!NOTE]
> Pip installation might encounter issues due to the dependencies mentioned above, especially on systems without `zlib.h` or on Windows.
> `sudo apt install gcc zlib1g zlib1g-dev` (Ubuntu)  
> `brew install gcc zlib` (macOS)


## Report other troubles

Please use [GitHub Issues](https://github.com/akikuno/DAJIN2/issues) for all reporting purposes.  
