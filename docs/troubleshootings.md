# Troubleshootings


## Installation

### Installation of `mappy` via pip was failed

These packages are needed to install `mappy`.

```bash
sudo apt install gcc build-essential python3-dev
pip install --upgrade pip mappy
```

### Installation of `pysam` was failed

Python 3.11.0 fails to install `pysam` (at 2022-11-03).

Downgrade Python to 3.10 and install `pysam`

```bash
conda install python=3.10 pysam
```
