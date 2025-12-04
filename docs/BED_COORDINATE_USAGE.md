# BED File Coordinate Support

This document describes the new `--genome-coordinate` option added to DAJIN2.

## Overview

DAJIN2 now supports specifying genomic coordinates via BED files using the `-b/--bed` option (or `--genome-coordinate` for backward compatibility). This allows users to provide precise genomic coordinates without needing to rely on UCSC BLAT lookup, and enables offline operation when the external servers (UCSC Genome Browser and GGGENOME) are unavailable.

## Usage

### Command Line Interface

```bash
DAJIN2 --sample sample/ --control control/ --allele alleles.fa --name experiment1 --bed coordinates.bed
```

### Priority Rules

1. If `-b/--bed` is specified, it takes **precedence** over `--genome`
2. If both are specified, the BED file coordinates are used, but the genome ID from `--genome` is included in the metadata
3. If only `--genome` is specified, UCSC lookup is used as before
4. If neither is specified, analysis proceeds without genomic coordinates

### BED File Format

DAJIN2 accepts BED6 format files with the following requirements:

```
chr1	1000000	1001000	mm39	248956422	+
chr2	2000000	2001000	mm39	242193529	-
```

#### Column Definitions
- **Column 1**: Chromosome name (e.g., `chr1`, `chr2`, `chrX`, or `1`, `2`, `X`)
- **Column 2**: Start position (0-based, inclusive)
- **Column 3**: End position (0-based, exclusive)
- **Column 4**: Name (**genome ID**, e.g., `mm39`, `hg38`)
- **Column 5**: Score (**chromosome size for proper IGV visualization**, e.g., `248956422` for chr1)
- **Column 6**: Strand (**REQUIRED**, `+` or `-`, **must match FASTA allele orientation**)

> [!NOTE]  
> For the score field (column 5), please enter the size of the chromosome specified in column 1.  
> While the original BED format limits scores to 1000, DAJIN2 accepts **chromosome sizes without any issue**.

> [!IMPORTANT]  
> **Strand orientation must match**. The strand field (column 6: `+` or `-`) in your BED file **must match the strand orientation of your FASTA allele sequences**.  
> - If your FASTA allele sequence is on the **forward strand** (5' to 3'), use `+` in the BED file  
> - If your FASTA allele sequence is on the **reverse strand** (3' to 5'), use `-` in the BED file

**Examples:**
```
chr1    1000000    1001000    mm39    248956422    +
chr2    2000000    2001000    mm39    242193529    -
chrX    5000000    5001000    mm39    156040895    +
```

**How to find chromosome sizes:**
- Human (hg38): https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&chromInfoPage=
- Mouse (mm39): https://genome.ucsc.edu/cgi-bin/hgTracks?db=mm39&chromInfoPage=
- Or use: `samtools faidx genome.fa` to get chromosome sizes from FASTA index

#### Strand Handling
- **Strand field (column 6)** is **REQUIRED** for proper sequence orientation
- **`+` strand**: Sequence is on the forward strand (5' to 3')
- **`-` strand**: Sequence is on the reverse strand, will be reverse complemented
- **Required format**: DAJIN2 requires BED6 format with strand information (`+` or `-`)
- **IGV visualization**: Strand information is preserved in BAM files for genome browser display
- **Sequence processing**: DAJIN2 automatically applies reverse complement for minus strand regions

### File Extensions

Supported file extensions:
- `.bed`
- `.bed.gz`

### Multi-Interval BED Files

- If the BED file contains multiple intervals, DAJIN2 uses the **first interval** for analysis
- A warning is logged when multiple intervals are detected
- Future versions may support multi-region analysis

## Batch Mode Support

The `-b/--bed` option is also supported in batch mode:

### Batch File Format

```csv
sample,control,allele,name,genome,genome_coordinate
sample1/,control/,alleles.fa,exp1,hg38,coords1.bed
sample2/,control/,alleles.fa,exp2,,coords2.bed
sample3/,control/,alleles.fa,exp3,mm39,
```

### Batch Processing Rules

- Each row can specify its own BED file in the `genome_coordinate` column
- BED files are validated during batch file processing
- Same priority rules apply: `genome_coordinate` takes precedence over `genome`

## Examples

### Example 1: Basic BED File Usage

**coordinates.bed:**
```
chr1	140453136	140753136	mm39	248956422	+
```

**Command:**
```bash
DAJIN2 --sample test_sample/ --control test_control/ --allele test_alleles.fa --name test_analysis --bed coordinates.bed
```

### Example 2: Batch Mode with Mixed Coordinate Sources

**batch.csv:**
```csv
sample,control,allele,name,genome,bed
exp1/,ctrl/,alleles.fa,experiment1,hg38,
exp2/,ctrl/,alleles.fa,experiment2,,region2.bed
```

- `experiment1`: Uses `hg38`
- `experiment2`: Uses `region2.bed` without genome metadata

## Technical Details

### Coordinate System

- **Input**: BED format uses 0-based start, 1-based end (half-open interval)
- **Internal**: DAJIN2 maintains 0-based, half-open intervals `[start, end)`
- **Output**: Mutation coordinates in CSV use 0-based positions

### Validation

BED files are validated for:
- File existence and proper extension (`.bed` or `.bed.gz`)
- BED6 format with strand information (6 columns required)
- Valid coordinate ranges (start < end, start >= 0)
- Proper data types (integer coordinates)
- Chromosome size (4th column must be a positive integer)
- Strand information (required, must be `+` or `-` in field 6)

### Error Handling

Common errors and solutions:

**Error**: `BED file must have .bed or .bed.gz extension`
- **Solution**: Rename your file to have the correct extension

**Error**: `DAJIN2 requires BED6 format with strand information`
- **Solution**: Ensure your BED file has 6 columns including strand (+/-) in the 6th column

**Error**: `Invalid chromosome size format at line X: 'Y' (must be an integer)`
- **Solution**: Set the 4th column to a positive integer representing chromosome size (e.g., 248956422 for chr1)

**Error**: `Invalid chromosome size at line X: Y (must be a positive integer)`
- **Solution**: Ensure the 4th column contains a positive integer greater than 0

**Error**: `Invalid or missing strand at line X: 'Y' (must be '+' or '-')`
- **Solution**: Set the 6th column to either `+` or `-` for strand information

**Error**: `Invalid end position at line X: Y (must be > start)`
- **Solution**: Check that end position is greater than start position

**Error**: `No valid intervals found in BED file`
- **Solution**: Ensure the file contains at least one valid interval (not just comments)

## Integration with Existing Features

### Mutation Reporting

When genomic coordinates are provided via BED file:
- Mutation CSV includes genome assembly information
- Start/End columns contain absolute genomic positions
- HTML visualization displays genomic context

### IGV Integration

- Genomic coordinates enable proper IGV.js visualization
- Mutations are displayed in genomic context
- Navigation uses genome coordinates

## Backward Compatibility

- All existing `--genome` functionality is preserved
- Existing batch files without `genome_coordinate` column continue to work
- No changes to output formats or file structures

