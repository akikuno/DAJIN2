# BED File Coordinate Support - Issue #26

This document describes the new `--genome-coordinate` option added to DAJIN2.

## Overview

DAJIN2 now supports specifying genomic coordinates via BED files using the `--genome-coordinate` option. This allows users to provide precise genomic coordinates without needing to rely on UCSC BLAT lookup.

## Usage

### Command Line Interface

```bash
# Using BED file only
DAJIN2 --sample sample/ --control control/ --allele alleles.fa --name experiment1 --genome-coordinate coordinates.bed

# Using BED file with genome ID (BED takes precedence)
DAJIN2 --sample sample/ --control control/ --allele alleles.fa --name experiment1 --genome hg38 --genome-coordinate coordinates.bed
```

### Priority Rules

1. If `--genome-coordinate` is specified, it takes **precedence** over `--genome`
2. If both are specified, the BED file coordinates are used, but the genome ID from `--genome` is included in the metadata
3. If only `--genome` is specified, UCSC lookup is used as before
4. If neither is specified, analysis proceeds without genomic coordinates

### BED File Format

DAJIN2 accepts standard BED format files with the following requirements:

#### Minimum Format (3 columns)
```
chr1	1000000	1001000
chr2	2000000	2001000
```

#### Extended Format (6 columns)
```
chr1	1000000	1001000	248956422	0	+
chr2	2000000	2001000	242193529	0	-
```

#### Column Definitions
- **Column 1**: Chromosome name (e.g., `chr1`, `chr2`, `chrX`)
- **Column 2**: Start position (0-based, inclusive)
- **Column 3**: End position (0-based, exclusive)
- **Column 4**: Feature name (**RECOMMENDED**: use chromosome size, e.g., `248956422` for chr1)
- **Column 5**: Score (optional, typically 0)
- **Column 6**: Strand (**REQUIRED**, `+` or `-`)

#### Chromosome Size Recommendation
For optimal IGV visualization and proper coordinate handling, it is **strongly recommended** to specify the chromosome size in the 4th column (feature name field). This helps DAJIN2 properly handle genomic coordinates and improves visualization in genome browsers.

**Examples:**
```
chr1    1000000    1001000    248956422    0    +
chr2    2000000    2001000    242193529    0    -
chrX    5000000    5001000    156040895    0    +
```

**How to find chromosome sizes:**
- Human (hg38): https://genome.ucsc.edu/cgi-bin/hgTracks?chromInfoPage=
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
chr7	140453136	140753136	CFTR_locus	0	+
```

**Command:**
```bash
DAJIN2 --sample cf_sample/ --control cf_control/ --allele cftr_alleles.fa --name cftr_analysis --genome-coordinate coordinates.bed
```

### Example 2: BED File with Genome ID

**Command:**
```bash
DAJIN2 --sample sample/ --control control/ --allele alleles.fa --name analysis --genome hg38 --bed target.bed
```

This command:
1. Uses coordinates from `target.bed`
2. Includes `hg38` in the metadata
3. Logs: "Using BED file coordinates: target.bed"

### Example 3: Batch Mode with Mixed Coordinate Sources

**batch.csv:**
```csv
sample,control,allele,name,genome,bed
exp1/,ctrl/,alleles.fa,experiment1,hg38,region1.bed
exp2/,ctrl/,alleles.fa,experiment2,hg38,
exp3/,ctrl/,alleles.fa,experiment3,,region3.bed
```

- `experiment1`: Uses `region1.bed` with `hg38` metadata
- `experiment2`: Uses UCSC lookup for `hg38`
- `experiment3`: Uses `region3.bed` without genome metadata

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

## Testing

The BED coordinate feature includes comprehensive test coverage:

```bash
# Run BED handler tests
pytest tests/src/utils/test_bed_handler.py -v

# Run input validator tests
pytest tests/src/utils/test_input_validator.py::TestValidateBedFileAndGetCoordinates -v
```

## Future Enhancements

Potential future improvements:
- Support for multiple genomic regions in a single analysis
- BED file output format for detected mutations
- Integration with additional genome browsers
- Support for BED detail formats (BED12, etc.)