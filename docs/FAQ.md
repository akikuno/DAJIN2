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

## Why is the read count of the Control sample lower in the output BAM file?

In DAJIN2 (version 0.5.4 and later), if the read count of the Control sample exceeds 10,000, it is randomly subsampled to 10,000 reads to reduce computational load. As a result, the BAM file contains fewer reads than the original count due to this filtering.  
It has been confirmed that limiting the read count to 10,000 does not affect the analysis results.
