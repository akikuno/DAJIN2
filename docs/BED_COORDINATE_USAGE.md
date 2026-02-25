# Specifying Genome Coordinates with BED Files

This document explains the BED file required by the `-b/--bed` option.

---

## Overview

DAJIN2 currently supports genome coordinate specification via BED files using the `-b/--bed` option (`--genome-coordinate` is also available for backward compatibility). This allows DAJIN2 to run offline even when external servers (UCSC Genome Browser and GGGENOME) are unavailable.

> [!IMPORTANT]
> The `--genome` option is still valid, but we recommend `-b/--bed` because communication errors with external servers may occur.

### Priority Rules

1. If `-b/--bed` is specified, it takes priority over `--genome`
2. If both are specified, coordinates from the BED file are used
3. If only `--genome` is specified, genome coordinates are obtained through communication with external servers (UCSC Genome Browser and GGGENOME), as before
4. If neither is specified, analysis runs without genome coordinates

### BED File Format

DAJIN2 accepts BED6-format files that meet the following requirements.

````
chr1	1000000	1001000	mm39	195154279	+
````

#### Column Definitions

- **Column 1**: Chromosome name (for example, `chr1`, `chr2`, `chrX`, or `1`, `2`, `X`)
- **Column 2**: Start position (0-based index)
- **Column 3**: End position (0-based index)
- **Column 4**: Name (**genome ID**, for example, `mm39`, `hg38`)
- **Column 5**: Score (**chromosome size for proper IGV visualization**, for example, `248956422` for chr1)
- **Column 6**: Strand (`+` or `-`, **must match FASTA allele orientation**)

> [!NOTE]
> In the score field (Column 5), enter the chromosome size of the chromosome specified in Column 1.
> In the original BED format, score is limited to 1000, but DAJIN2 accepts **chromosome sizes without issue**.

> [!IMPORTANT]
> **Strand direction must match**. The strand field in the BED file (Column 6: `+` or `-`) **must match the strand direction of the FASTA allele sequence**.
> - If the FASTA allele sequence is on the **forward strand** (5'→3'), use `+` in the BED file
> - If the FASTA allele sequence is on the **reverse strand** (3'→5'), use `-` in the BED file

**How to check chromosome size:**

- Human (hg38): https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&chromInfoPage=
- Mouse (mm39): https://genome.ucsc.edu/cgi-bin/hgTracks?db=mm39&chromInfoPage=
- Or use `samtools faidx genome.fa` to get chromosome sizes from a FASTA index

#### Strand Handling

- The **strand field (Column 6)** is required for proper sequence orientation
- **`+` strand**: Sequence is on the forward strand (5'→3')
- **`-` strand**: Sequence is on the reverse strand and converted to its reverse complement
- **Required format**: DAJIN2 requires BED6 format including strand information (`+` or `-`)
- **IGV visualization**: Strand information is preserved in BAM files and used for genome browser display
- **Sequence processing**: DAJIN2 automatically applies reverse complement to minus-strand regions

### File Extensions

Supported file extensions:

- `.bed`
- `.bed.gz`

### BED Files Containing Multiple Intervals

- If a BED file contains multiple intervals, DAJIN2 uses the **first interval** for analysis
- A warning log is output when multiple intervals are detected
- Future versions may support multi-region analysis

---

## Usage

### Single-Sample Analysis

````bash
DAJIN2 --sample sample/ \
    --control control/ \
    --allele alleles.fa \
    --name experiment1 \
    --bed coordinates.bed
````

### Batch Mode Support

The `-b/--bed` option is also supported in batch mode.

### Batch File Format

````csv
sample,control,allele,name,genome,bed
sample1/,control/,alleles.fa,exp1,hg38,exp1.bed
sample2/,control/,alleles.fa,exp2,,exp2.bed
````

- You can specify a custom BED file in the `bed` column for each row
- The same priority rule applies: `bed` takes precedence over `genome`

---

## Error Handling

Common errors and solutions:

**Error**: `BED file must have .bed or .bed.gz extension`

- **Solution**: Rename the file to use a valid extension

**Error**: `DAJIN2 requires BED6 format with strand information`

- **Solution**: Ensure the BED file has 6 columns and that Column 6 contains strand information (`+/-`)

**Error**: `Invalid chromosome size format at line X: 'Y' (must be an integer)`

- **Solution**: Set Column 5 to a positive integer representing chromosome size

**Error**: `Invalid chromosome size at line X: Y (must be a positive integer)`

- **Solution**: Ensure Column 5 contains a positive integer greater than 0

**Error**: `Invalid or missing strand at line X: 'Y' (must be '+' or '-')`

- **Solution**: Set Column 6 to either `+` or `-` as strand information

**Error**: `Invalid end position at line X: Y (must be > start)`

- **Solution**: Ensure the end position is greater than the start position

**Error**: `No valid intervals found in BED file`

- **Solution**: Ensure the BED file contains at least one valid genome coordinate interval

---

## Backward Compatibility

- All existing `--genome` functionality is preserved
- Existing batch files without a `genome_coordinate` column continue to work
- No changes are made to output format or file structure
