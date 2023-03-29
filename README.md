# DAJIN2 is now under development👷

I aim to release a stable version of DAJIN2 by ~~March~~ August 2023:crossed_fingers:

## New planned features


### Experiments

+ [x] Nanopore R9.4
+ [x] Nanopore R10.4
+ [ ] PacBio
+ [ ] nCATs
+ [ ] base-editing

### Tested Design

+ [x] point mutation
+ [x] 2-cut knockout
+ [x] flox knockin
+ [ ] 1-cut knockout
+ [ ] large insertion knockin
+ [ ] humanized exon knockin
+ [ ] multi-allele phasing

### Functions

+ [ ] output VCF
+ [ ] multi-processing
+ [x] error correction
  - repetitive deletion
  - knockin sequence
+ [x] trim microhomology in BAM

### Usability

+ [ ] documentation
+ [ ] easy installation
  + [x] pypi
  + [ ] bioconda
  + [ ] docker image
+ [ ] support low-/middle-spec laptop
+ [x] support macOS
+ [x] support samples without a reference genome
+ [x] support visualization using `igv.js` via `DAJIN2 view` command
+ [x] support graphical user interface (GUI) using `flask` via `DAJIN2 gui` command

### Maintainability

+ [ ] tests
+ [x] minimize dependencies
  + [x] Python, R, Bash -> Python
  + [x] samtools -> pysam
  + [x] minimap2 -> mappy
  + [x] remove NanoSim
  + [x] remove Tensorflow (GPU computation)

