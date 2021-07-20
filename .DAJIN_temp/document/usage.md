Usage:

  DAJIN [options] -a <alleles.fasta> -c <control.fastq> -s <sample.fastq>

Example:

  DAJIN \
    -a example/alleles.fa \
    -c example/fastq/barcode01.fq.gz \
    -s example/fastq/barcode02.fq.gz \
    -g mm10 \
    -o DAJIN_example \
    -t 4

Mandatory arguments:

  -a|--alleles <string>
    Name of a multi-fasta file describing DNA sequence of each alleles in specific regions of interest.

  -c|--control <string>
    Name of a fastq file from a control sample.
    Acceptable formats include gzipped FASTQ and FASTQ.

  -s|--sample <string>
    Name of a fastq file from an experimental sample.
    Acceptable formats include gzipped FASTQ and FASTQ.

Optional arguments (defaults in parentheses):

  -g|--genome <string>
    Name of UCSC genome releases.
    The available genomes are listed at https://genome.ucsc.edu/FAQ/FAQreleases.html#release1

  -t|--threads <integer>
    Number of threads to use (1)

  -o|--output <string>
    Output directory (DAJIN_results)

See the further information: <https://github.com/akikuno/DAJIN2>

