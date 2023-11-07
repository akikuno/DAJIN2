# Troubleshooting

## Installation Issues

### Failure in Installing `mappy` with pip

To successfully install `mappy`, the following packages are required:

```bash
sudo apt install gcc build-essential python3-dev
pip install --upgrade pip
pip install mappy
```

### Failure in Installing `pysam` with Python 3.11.0

As of 2023-11-08, Python 3.11.0 cannot install `pysam`.  
It's recommended to downgrade Python to version 3.10 and then attempt the installation:

```bash
conda install python=3.10 pysam --yes
```

### Failure in Installing `weasyprint` on macOS

On macOS, `weasyprint` should be installed using `homebrew`:

```bash
brew install weasyprint
```
