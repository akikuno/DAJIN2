# Troubleshooting Guide

## Installation

### Prerequisites

Before installing `DAJIN2`, please ensure that your system meets the following requirements:

- Python >=3.7
- Unix-like environment (Linux, macOS, WSL, etc.)
- [conda](https://docs.conda.io/en/latest/) or [mamba](https://mamba.readthedocs.io/en/latest/) is highly recommended for managing dependencies
- If using pip, access to administrative privileges or the ability to install packages globally

### Dependencies

`DAJIN2` depends on the `pysam` package, which in turn requires `htslib`. These dependencies are critical for the functionality of `DAJIN2` but can pose installation challenges:

- `htslib` requires `zlib.h`, which may not be available on all systems by default.
- As of the latest update, `htslib v1.18`, there is no support for Windows via Bioconda.

### Recommended Installation Method

#### Using Conda/Mamba

We strongly recommend using Conda or Mamba for installation, as they efficiently manage the complex dependencies of `pysam` and `htslib`:

1. Install Conda or Mamba if you haven't already. You can download Conda from [here](https://docs.conda.io/en/latest/miniconda.html).

2. Create a new environment (optional but recommended):
```bash
conda create -n DAJIN2-env
conda activate DAJIN2-env
```

3. Install `DAJIN2`:
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


## Report other troubles

Please use [GitHub Issues](https://github.com/akikuno/DAJIN2/issues) for all reporting purposes.  
