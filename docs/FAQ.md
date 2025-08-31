# Frequently Asked Questions

## How many reads are necessary?

**We recommend at least 1,000 reads.**  
With 1,000 reads, it is possible to detect point mutation alleles with a frequency of 1%, ensuring high-precision analysis. However, if the target is an indel of more than a few tens of bases, or if the expected allele frequency is 5% or higher, detection is possible with fewer reads (~500 reads).

## What is the recommended read length for analysis?

**We recommend lengths below 10kb.**  
If the length is below 10kb, when reading PCR amplicons with Nanopore, the reads will uniformly cover the target region. It is possible to obtain PCR amplicons up to approximately 15kb, but this may result in uneven coverage of the target region, potentially reducing analysis accuracy.

## Can data from platforms other than Nanopore (e.g., PacBio or NGS) be analyzed?

**Yes, it is possible.**  
DAJIN2 accepts common file formats (FASTA, FASTQ, BAM) as input, allowing the analysis of data from platforms other than Nanopore. However, since we do not have experience using DAJIN2 with non-Nanopore data, please contact us [here](https://github.com/akikuno/DAJIN2/issues/new/choose) if you encounter any issues.

## Some samples are not processed when `--threads` is set to 1 or higher in batch mode

**Please update DAJIN2 to the latest version to fix the bug.**  

```bash
conda update DAJIN2 -y
```

or

```bash
pip install -U DAJIN2
```

## Why is the read count of the control sample low in the output bam file?

In DAJIN2 (version 0.5.4 and later), if the read count of the Control sample exceeds 100,000, the number of reads is randomly limited to 100,000 to reduce computational load.  
As a result, the BAM file contains filtered reads, leading to a lower read count than the original.  
It has been confirmed that limiting the read count to 100,000 does not affect the analysis results.  

## What should I do if I get an error about GGGENOME or UCSC server being unavailable?

When using the `-g/--genome` option, DAJIN2 attempts to access the GGGENOME and UCSC Goldenpath servers to retrieve genomic coordinates.  
If you encounter error messages indicating that these servers are down or unavailable, consider using the `-b/--bed` option instead.  
The BED file option allows you to specify genomic coordinates directly without requiring external server access, enabling offline operation.  

For details on how to create and use BED files, see [README.md](https://github.com/akikuno/DAJIN2#using-bed-files-for-genomic-coordinates).  

## How can I use reference genomes other than those provided by UCSC?

You can specify any reference genome using the `-b/--bed` option.

For details on how to create and use BED files, see [README.md](https://github.com/akikuno/DAJIN2#using-bed-files-for-genomic-coordinates).  
