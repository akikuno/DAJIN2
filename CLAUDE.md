# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

DAJIN2 is a comprehensive genotyping tool for analyzing genome-edited samples using nanopore target sequencing. It detects mutations from point mutations to structural variants, classifies mosaic alleles, and provides intuitive HTML visualizations of results.

## Development Commands

### Installation and Setup
```bash
# Development installation
pip install -e .

# Install with conda environment
conda create -y -n env-dajin2 python=3.12 pip pytest
conda activate env-dajin2
pip install -e .
```

### Testing
```bash
# Run all tests
pytest tests/ -p no:warnings -vv

# Run tests with slow tests included
pytest tests/ --runslow -p no:warnings -vv

# Run specific test module
pytest tests/src/preprocess/test_sv_detector.py -v

# Run single test function
pytest tests/src/utils/test_fastx_handler.py::test_specific_function -v
```

### Code Quality
```bash
# Format and lint code (using ruff configuration from pyproject.toml)
ruff format .
ruff check .

# Type checking
mypy src/DAJIN2/
```

### Running DAJIN2
```bash
# Single sample analysis
DAJIN2 --sample path/to/sample --control path/to/control --allele sequences.fa --name output_name

# Batch processing
DAJIN2 batch --file batch.csv --threads 4

# GUI mode
DAJIN2 gui

# View results with IGV.js
DAJIN2 view --name output_name
```

## Architecture Overview

### Core Processing Pipeline
The main processing flow follows this pattern:
1. **Input Validation & Formatting** (`preprocess/input_formatter.py`)
2. **Control Processing** (`core.execute_control()`)
3. **Sample Processing** (`core.execute_sample()`)
4. **Report Generation** (`report_generator.py`)

### Major Components

#### Core Processing (`src/DAJIN2/core/`)
- **preprocess/**: Input handling, mapping, SV detection, mutation extraction
- **clustering/**: K-mer based clustering, label management, strand bias handling
- **consensus/**: Consensus sequence generation, SV annotation
- **classification/**: Allele classification and merging
- **report/**: BAM, FASTA, HTML, and mutation info export

#### Key Processing Steps
1. **Preprocessing Pipeline**:
   - `mapping.py`: Sequence alignment using mappy
   - `sv_detector.py`: Structural variant detection
   - `mutation_extractor.py`: Extract mutations from alignments
   - `midsv_caller.py`: Mid-SV calling using external midsv library

2. **Clustering Pipeline**:
   - `kmer_generator.py`: Generate k-mers for similarity analysis
   - `clustering.py`: Cluster sequences by similarity
   - `label_updator.py`: Update cluster labels based on mutations

3. **Consensus Pipeline**:
   - `consensus.py`: Generate consensus sequences for each cluster
   - `sv_annotator.py`: Annotate structural variants in consensus

#### Utilities (`src/DAJIN2/utils/`)
- `fastx_handler.py`: FASTQ/FASTA file operations
- `sam_handler.py`: SAM/BAM file operations using pysam
- `cssplits_handler.py`: Handle CS tags for mutations
- `multiprocess.py`: Parallel processing utilities

### Data Flow Architecture
The system uses a temporary directory structure under `DAJIN_Results/.tempdir/` during processing:
- Input files are preprocessed and stored in standardized formats
- Intermediate results are cached to enable resuming interrupted runs
- Final results are exported to `DAJIN_Results/{output_name}/`

### Configuration
- Global configuration is managed through `utils/config.py`
- Logging is configured per-run with timestamps
- Thread management is handled automatically based on available cores

## Testing Architecture

Tests are organized to mirror the source code structure:
- `tests/src/` contains unit tests for each module
- `tests/data/` contains test datasets for integration testing
- `conftest.py` configures pytest with slow test markers
- Use `@pytest.mark.slow` for tests that take significant time

## Important Implementation Details

### Memory Management
- Large datasets are processed in chunks to manage memory usage
- Temporary files are automatically cleaned up unless `--debug` flag is used
- Environment variables control threading to prevent oversubscription

### External Dependencies
- **mappy**: For sequence alignment
- **pysam**: For BAM file operations
- **midsv**: Custom library for structural variant calling
- **cstag**: For handling CS tags in SAM files

### Multi-threading Architecture
- Batch mode processes multiple samples in parallel
- Individual sample processing uses configurable thread counts
- Thread counts are validated and adjusted based on system capabilities

### Input/Output Formats
- Supports FASTQ/FASTA (compressed/uncompressed) and BAM input
- Outputs include BAM files, FASTA sequences, HTML visualizations, CSV mutation tables
- HTML reports use embedded IGV.js for interactive visualization

## Version and Dependencies

Current version: 0.6.2
- Python 3.9-3.12 required
- Uses Poetry for dependency management
- All dependencies specified in pyproject.toml
- Supports Unix-based environments (Linux, macOS, WSL2)
