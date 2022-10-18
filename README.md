# DAJIN2 is now under development👷

I aim to release a stable version of DAJIN2 by March 2023:crossed_fingers:
## New planned features


### Wet experiments

- [ ] support nCATs
- [ ] support PacBio

### Functions

- [ ] output VCF
- [ ] multi-allele phasing
- [x] trim microhomology in BAM

### Usability

- [ ] easy installation
  - [ ] PyPI
  - [ ] BioConda
  - [ ] docker image
- [ ] support low-/middle-spec laptop
- [x] support macOS
- [x] support samples without a reference genome
- [x] support visualization using `igv.js` via `DAJIN2 view` command
- [x] support graphical user interface using `flask` via `DAJIN2 gui` command

### Maintainability

- [ ] tests
- [x] minimize dependencies
  - [x] Python, R, Bash -> Python
  - [x] samtools -> pysam
  - [x] minimap2 -> mappy
  - [x] remove NanoSim
  - [x] remove Tensorflow (GPU computation)

