<!-- TEMPLATE
# v0.0.0 (yyyy-mm-dd)
## 💥 Breaking
## 📝 Documentation
## 🚀 Performance
## 🌟 New Features
## 🐛 Bug Fixes
## 🔧 Maintenance
## ⛔️ Deprecated
+ commitMessage. Issue [#XX](https://github.com/akikuno/DAJIN2/issues/XX) [[Commit Detail](https://github.com/akikuno/DAJIN2/commit/xxxxx)]
-->

<!-- ############################################################# # -->

<!-- ############################################################# # -->


## 🌟 New Features

+ Added `--no-filter` option to detect rare mutations. See Issue: [#83](https://github.com/akikuno/DAJIN2/issues/83)

+ Added `-b/--bed` option to specify a BED file when using genomes other than UCSC reference genomes. See Issue: [#26](https://github.com/akikuno/DAJIN2/issues/26)

+ Added a feature to display an explicit error message with a UCSC genome browser link when the input FASTA sequence is not found in the reference genome. See Issue: [#26](https://github.com/akikuno/DAJIN2/issues/26)

+ Removed the `ValueError` that occurred when identical sequences were found. This is to support short-read data, where such cases are expected. [[Commit Detail](https://github.com/akikuno/DAJIN2/commit/fd7518357779f2f2ebea1dd6a43c7ddfe0b4c5a3#diff-4b9f2a326855933258d70bf13c107eed755b5140240b5b30bc9ca05de397ceeb)]

## 🔧 Maintenance

+ Updated Python support from 3.9 to 3.12 due to dependencies in `pysam` and `mappy`. See Issue: [#101](https://github.com/akikuno/DAJIN2/issues/101)


-------------------------------------------------------------

# Past Releases

<!-- <details>
<summary> v0.X.X (2025-MM-DD) </summary>

</details> -->

<details>

<summary> v0.6.2 (2025-06-07) </summary>

## 📝 Documentation

+ Update README to clarify hardware requirements Issue: [#91](https://github.com/akikuno/DAJIN2/issues/91), [#92](https://github.com/akikuno/DAJIN2/issues/92) [[Commit Detail](https://github.com/akikuno/DAJIN2/commit/a7995727dc35944d51c318e6ae0e79ded992c684)]


## 🐛 Bug Fixes

+ Fix CSV header validation error when processing files with BOM (Byte Order Mark). Issue [#88](https://github.com/akikuno/DAJIN2/issues/88). [[Commit Detail](https://github.com/akikuno/DAJIN2/commit/2bbf405a53d1fe431999cd368f55a3f0e422622b)]

+ The strand bias detection was too lenient, so the following revisions were made. Issue [#89](https://github.com/akikuno/DAJIN2/issues/89):
  + In the `preprocessing` step, strand bias was calculated at the nucleotide level, and nucleotides showing bias were excluded. [[Commit Detail](https://github.com/akikuno/DAJIN2/commit/c69dba52c023731aeb6d5f3c75b7550fcc426961)]
  + In the `clustering` step, strand bias was calculated at the cluster level, and clusters showing bias were excluded. [[Commit Detail](https://github.com/akikuno/DAJIN2/commit/cb1286250cd04ab592d81bf532558f29f211f91f)]

+ The order of arguments in `fastx_handler.convert_bam_to_fastq` was reversed, so it has been corrected. Issue [#94](https://github.com/akikuno/DAJIN2/issues/94) [[Commit Detail](https://github.com/akikuno/DAJIN2/commit/9498272d4afad04df5e2b662b4fedd8c1371024f)]


## 🔧 Maintenance

+ Shorten DAJIN2 log filename by replacing long UUID with microseconds for visual clarity. Issue [#95](https://github.com/akikuno/DAJIN2/issues/95) [[Commit Detail](https://github.com/akikuno/DAJIN2/commit/aceff9aec1f0146434f14597df2f77c508b0726f)]


</details>


<!-- ############################################################# # -->

<details>
<summary> v0.6.1 (2025-03-18) </summary>

## 🚀 Performance

- Use `BisectingKMeans` instead of `AgglomerativeClustering` because `BisectingKMeans` can take a `spmatrix` as input, significantly reducing memory consumption. [[Commit Detail](https://github.com/akikuno/DAJIN2/commit/19fe3549584ee1b1c8ccb67c3d364434f5ad392c)]

## 📝 Documentation
+ Specify the Range of Bases to Be Recorded in the FASTA File. Issue #78 [[Commit Detail](https://github.com/akikuno/DAJIN2/commit/a8955c5dc38f28279a708b485466081c5f39aa3e)]

## 🔧 Maintenance

+ Explicitly unify the line endings of text files in DAJIN_Reports to `LF`. [[Commit Detail](https://github.com/akikuno/DAJIN2/commit/695d67a7b6de3f29381bd4299fddc6106028d5c2)]

+ Upgrade to `pandas = ">=2.0.0"` because the argument specification for line terminator was changed to `lineterminator` in pandas >=1.5. [[Commit Detail](https://github.com/akikuno/DAJIN2/commit/65b421569e8be93e06c496827c8de1c62b9f47e0)]

+ Sort MUTATION_INFO by Allele ID. Issue [#79](https://github.com/akikuno/DAJIN2/issues/79) [[Commit Detail](https://github.com/akikuno/DAJIN2/commit/891096d963d408e454d26f3eab26a5cee6426b3a)]

+ Add `sv_annotator` to reflect SV alleles in consensus midsv tags. [[Commit Detail](https://github.com/akikuno/DAJIN2/commit/1f8a89d50c53427ae4166e5d48f57466441c2f50)]

+ Refactoring `annotate_insertions_within_deletion`: Previously, a similar function existed in `cssplit_handler`, but since this function is only called once during consensus, it has been moved to a dedicated module, `consensus.sv_annotator`. At the same time, the function has been simplified. [[Commit Detail](https://github.com/akikuno/DAJIN2/commit/88988db9cfeb44fd95f377a38eb434289ac812b5)]
 
## 🐛 Bug Fixes

+ Reflect Inversion Alleles When Flanked by Deletions at HTML. Issue #82 [[Commit Detail](https://github.com/akikuno/DAJIN2/commit/c7bb3cfacf724dd0a4298c432bdda878f5a72de4)]

+ Fix the issue where the SV length was reflected one base longer in deletion/inversion SV alleles. [[Commit Detail](https://github.com/akikuno/DAJIN2/commit/cb39c83d3bea0f3d507dbc4d22c5a58f6d3cee3f)]

+ Fix a bug where the silhouette score could not be calculated and resulted in an error when the sample and control were completely separated at a 1:1 ratio. [[Commit Detail](https://github.com/akikuno/DAJIN2/commit/c4f36483a03584f52299e2cba81d9eaf9d19425d)]

+ Correct the mislabeling of Deletion Allele as Insertion Allele. [[Commit Detail](https://github.com/akikuno/DAJIN2/commit/b9e38830e588823fc687363dee343f70fd5b6fae)]

+ Return the region containing the insertion sequence as a deletion sequence if the region flanked by deletions is determined to be an insertion sequence. Issue #86 [[Commit Detail](https://github.com/akikuno/DAJIN2/commit/eef0c0d1b2658cd22ebca1108495aed2adea6e20)]

+ Inversions are underlined since they can coexist with other mutations, while others are highlighted. Issue #84 [[Commit Detail](https://github.com/akikuno/DAJIN2/commit/ccdc98b97052dbbb27b8c75b3dd1b888620a4e11)]

+ Reflect the mutations (indel, substitution) within the inversion in HTML and MUTATION_LOCI. [[Commit Detail](https://github.com/akikuno/DAJIN2/commit/d0316c52cc129cf8d104cd1cce9c5eedff9898ce)]
+ 
</details>


<details>
<summary> v0.6.0 (2025-02-20) </summary>

## 💥 Breaking

+ Add `preprocess.sv_detector` to detect SV (Insertion/Deletion/Inversion) alleles. Issue #33 [[Commit Detail](https://github.com/akikuno/DAJIN2/commit/710e2a860abd9899d0b7378a58e1d672754c43db)]

+ Add `html_builder` to display SV alleles. Issue #31 [[Commit Detail](https://github.com/akikuno/DAJIN2/commit/1ba73f00ac72a5191749248469fab639c9e1e429)]

## 📝 Documentation

+ Upgrade Python version from 3.10 to 3.12 in README.md. Issue #74 [[Commit Detail](https://github.com/akikuno/DAJIN2/commit/9499666e78c669debe8e5f600dc6a8f82253bae8)]


## 🚀 Performance

+ Simplify feature extraction using `extract_n_features` to reduce computational costs. [[Commit Detail](https://github.com/akikuno/DAJIN2/commit/cc4665719f9eedfa0457cc07d1bc3ca2142c574e)]

+ To avoid overlooking minor alleles, the number of reads is increased from 10,000 to 100,000 during downsampling. [[Commit Detail](https://github.com/akikuno/DAJIN2/commit/7e3e08dec5d71769f917a58f63bf5af6df230ed4)]

## 🔧 Maintenance

+ Increase the SV allele number to at least two digits (e.g., `deletion01`). [[Commit Detail](https://github.com/akikuno/DAJIN2/commit/75367ea57fc1d3a4a9d86ea7e085a707da7ccdb3)]

+ Display the currently processing NAME in batch mode. [[Commit Detail](https://github.com/akikuno/DAJIN2/commit/c713f3198a9cc3f27854fec71c4c1da8959b4053)]

+ By appending a UUID to the log file, potential filename duplication can be prevented. [[Commit Detail](https://github.com/akikuno/DAJIN2/commit/41a8b183c47642ad2c2f171d02c779fd8f817c55)]

</details>

<details>
<summary> v0.5.6 (2024-12-06) </summary>

## 💥 Breaking

+ Support for PacBio HiFi reads. [[Commit Detail](https://github.com/akikuno/DAJIN2/commit/69773342d3bf157f28db013d62796be43ae297e7)]

+ Add `preprocess.sequence_error_handler` to exclude Nanopore sequence errors from the analysis. Issue: #60 
  + Initial commit [[Commit Detail](https://github.com/akikuno/DAJIN2/commit/329493dbfe6c0d03a6f8aadaab88911900f35dbb)]
  + Since most Nanopore sequencing errors occur due to read interruptions, `parse_midsv_from_csv` classifies entries as either Unknown or Other (M). [[Commit Detail](https://github.com/akikuno/DAJIN2/commit/809c22be361b0566e46ea94b5dec37b8a4659244)]
  + Instead of strategies like Cosine similarity or HDBSCAN, the Jaro-Winkler distance is explicitly used as a string similarity metric. Jaro-Winkler was chosen because Levenshtein would be too time-consuming. [[Commit Detail](https://github.com/akikuno/DAJIN2/commit/809c22be361b0566e46ea94b5dec37b8a4659244)]

+ Add `sr` presets to all execusions in `preprocess.mapping`. Issue: #55 [[Commit Detail](https://github.com/akikuno/DAJIN2/commit/682a20f3e71206fbd55369b5ff0dea799881aa67)]

+ Increase the sensitivity by lowering the mutation detection threshold from 0.5% to 0.1% to detect mutations around 0.75%. [[Commit Detail](https://github.com/akikuno/DAJIN2/commit/0e19752c1a4100d5a9121d54a563698642dc35c1)]

+ Use `AgglomerativeClustering` instead of Constrained KMeans because AgglomerativeClustering provides a more global clustering approach, and Constrained KMeans was not very useful due to the unreliability of its `min_cluster_size`. [[Commit Detail](https://github.com/akikuno/DAJIN2/commit/b81711e553edccaec4d6396cf940081163a18471)]

+ Output seqence error reads as `BAM/{name}/sequence_errors.bam`. Issue: #61 [[Commit Detail](https://github.com/akikuno/DAJIN2/commit/b2e717faa5ed17cdc18217d537cf49de6ca7c0b4)]

## 🚀 Performance

+ Downsampling the sample reads to a maximum of 10,000. Issue: #58 [[Commit Detail](https://github.com/akikuno/DAJIN2/commit/33c2120b59c80afe7f76165b2410b8ffe51410bd)]

## 🐛 Bug Fixes

+ Fix a bug where a element of dict with empty values was left behind after minor insertions were removed.  [[Commit Detail](https://github.com/akikuno/DAJIN2/commit/dcd66324999fbf2054b1310b879f81baf0fa7a92)]

## 🔧 Maintenance

+ With the end of security support for Python 3.8 in October 2024, we have updated DAJIN2 to support Python 3.9 or later. [[Commit Detail](https://github.com/akikuno/DAJIN2/commit/0967c463386a48639a614849b3d3e4453079c8b1)]

+ Replace typing.Generator to collections.abc.Iterator Since typing.Generator is deprecated. Issue: #53 [[Commit Detail](https://github.com/akikuno/DAJIN2/commit/f85964a4b8027b547b7b3e370b9e86ff8dda36be)]

+ Automatically retrieve version information using `importlib.metadata.version` Issue: #59 [[Commit Detail](https://github.com/akikuno/DAJIN2/commit/4cf75078b5e7f487b07650e934d63448bc3a328e)]

+ Move the FASTX IO processing to `utils.io`. Issue: #66 [[Commit Detail](https://github.com/akikuno/DAJIN2/commit/cca3db0b13ac53e082e2272f7ee7f593f905bd25)]

+ Add E2E tests in Github Actions. [[Commit Detail](https://github.com/akikuno/DAJIN2/commit/8fb93621ae7b9c4e867a68e9160c3295bfb0f872)]

</details>


<details>
<summary> v0.5.5 (2024-10-08) </summary>

## 📝 Documentation

+ Add `FAQ.md` and `FAQ_JP.md` to address the question: "Why is the read count of the Control sample lower in the output BAM file?". [[Commit Detail](https://github.com/akikuno/DAJIN2/commit/b238d21fbb7cd3330a147bdde65b726278447649)]

## 🔧 Maintenance

+ Integrating insertion and inversion detection: Issue #31
  + Add sv_handler [[Commit Detail](https://github.com/akikuno/DAJIN2/commit/d994d845b0b8ed0fa8affed7992f1d95bf163073)]

  + Modify arguments of `is_insertion` to `is_sv` [[Commit Detail](https://github.com/akikuno/DAJIN2/commit/f2d3dc4ca2dff60fc869fb1f5b6b08f54490b564)]

  + Remame `insertions_to_fasta.generate_insertions_fasta` to `insertion_detector.detect_insertions` because the function is not only for generating fasta files but also for generating csv tag. [[Commit Detail](https://github.com/akikuno/DAJIN2/commit/63c9d63bad627f529f272ea90c035e236f9dd1fb)]

+ Remove unused dependencies
  + `networkx`: Issue #49 [[Commit Detail](https://github.com/akikuno/DAJIN2/commit/524186bdce9e28d6357378d0baeb45670d2e22ed)]

</details>



<details>
<summary> v0.5.4 (2024-07-23) </summary>

## 💥 Breaking

+ Use simulated annealing to optimize cluster assignments in `clustering.constrained_kmenas` [[Commit Detail](https://github.com/akikuno/DAJIN2/commit/b07b626c1def93022e79840e1e6e393fa400cefb)]
  + Since `ortools` is not installable on osx-arm64 in Bioconda, I implemented alternative smethods to calculate min_cost_flow.

+ Change the criteria for terminating clustering. [[Commit Detail](https://github.com/akikuno/DAJIN2/commit/db6ec7245d0d1a7ff7204574cffdfd945ee5e854)]
  + The following termination criteria have been added:
    - Minimum cluster size is less than or equal to 0.5% of the sample's read number.
    - Decrease in the proportion of samples with a silhouette score of 0.25 or higher.
  + The following termination criterion has been removed:
    - Adjusted Rand Index >= 0.95, as it led to early termination when minor clusters were generated.

+ The threshold for `clustering.strand bias` determination has been loosened. [[Commit Detail](https://github.com/akikuno/DAJIN2/commit/5bbaa7d363bce03d6fbd4ba7fdf1c00e938d9809)]
  + This adjustment addresses cases like `+:13, -:2` (0.87) observed in `example_flox/flox-1nt-deletion`.
  + Since the minor allele is particularly susceptible, further adjustments may be necessary in the future.

## 🌟 New Features

+ Support for Apple Silicon (osx-arm64) in Bioconda. Issue: #46

</details>

<details>
<summary> v0.5.3 (2024-07-16) </summary>

## 💥 Breaking

- Update `clustering.clustering`: Use Constrained Kmeans clustering to address the issue of cluster imbalance where extremely minor clusters were preferentially separated. Set `min_cluster_size` to 0.5% of the sample read count. [[Commit Detail](https://github.com/akikuno/DAJIN2/commit/c1b14e73d8a95fdb39e510a7a90e501d596b7f3a)]
  - As a result, `clustering.label_merger.py` is no longer needed and has been removed.

- Update `consensus.call_consensus`: For mutations determined to be sequence errors, we previously replaced them with unknown (`N`), but this `N` had low interpretability. Therefore, mutations that DAJIN2 determines to be sequence errors will now be assigned the same base as the reference genome. [[Commit Detail](https://github.com/akikuno/DAJIN2/commit/1f46215ae7054c4da088c638ad82e41dd0dc7227)]

## 🐛 Bug Fixes

- Due to a bias in `classifiler.calc_match` where alleles with shorter sequences were prioritized, the operation of dividing by sequence length has been removed. [[Commit Detail](https://github.com/akikuno/DAJIN2/commit/fa6fbd5a7f9693df3b067a3041df42198a0d65b7)]

- Fix `preporcess.mapping.generate_sam` to perform alignments with `map-ont` and `splice` in addition to `sr` for sequence lengths of 500 bp or less, and select the optimal prefix from these alignments. Issue: #45 [[Commit Detail](https://github.com/akikuno/DAJIN2/commit/9e7fb93f3c7b74095d2afd08bf3fa0bc00e6f367)]
</details>


<details>
<summary> v0.5.2 (2024-07-08) </summary>

## 📝 Documentation

+ Add `FAQ.md` and `FAQ_JP.md` to provide answers to questions. [[Commit Detail](https://github.com/akikuno/DAJIN2/commit/c2217b006494ae73fda422a17edaf39fb97e8898)]

## 🌟 New Features

- Update `mutation_extractor` [[Commit Detail](https://github.com/akikuno/DAJIN2/commit/9444ee701ee52adeb6271552eff70667fb49b854)]
  - Simplified the logic of the `is_dissimilar_loci` if statement. Additionally, changed the threshold for determining a mutation in Consensus from 75% to 50% (to accommodate the insertion allele in Cas3 Tyr Barcode10).
  - Updated `detect_anomalies` to use MLPClassifier to detect mutations more flexibly and accurately compared to the previous threshold setting with MiniBatchKMeans.

## 🔧 Maintenance

+ Make DAJIN2 compatible with Python 3.11 and 3.12. Issue: #43 [[Commit Detail](https://github.com/akikuno/DAJIN2/commit/8da9118f5c0f584ed1ab12541d5e410d1b9f0da8)]
  + pysam and mappy builds with Python 3.11 and 3.12 are now available on Bioconda.

+ Update GitHub Actions to test with Python 3.11 and 3.12. Issue: #43 [[Commit Detail](https://github.com/akikuno/DAJIN2/commit/54df79e60b484da429c1cbf6f12b0c19196452cc)]

+ Resolve the B023 Function definition does not bind loop variable `alignment_lengths` issue. [[Commit Detail](https://github.com/akikuno/DAJIN2/commit/9c85d2f0410494a9b71d9905fad2f9e4efe30ed7)]

+ Add `question.yml` in GitHub Issue template. [[Commit Detail](https://github.com/akikuno/DAJIN2/commit/1172fddd34c382f92b6778d6f30fd733b458cc04)]


## 🐛 Bug Fixes

+ Update `cssplits_handler._get_index_of_large_deletions`: Modified to split large deletions when a match of 10 or more bases is found within the identified large deletion. Issue: #42 [[Commit Detail](https://github.com/akikuno/DAJIN2/commit/0c97a9b5fb8cad2ebdaf91b796eed3ce80f5eeee)]

</details>

<details>
<summary> v0.5.1 (2024-06-15) </summary>

## 🌟 New Features

+ Enable to accept additional file formats as an input. Issue: #37
  + FASTA [[Commit Detail](https://github.com/akikuno/DAJIN2/commit/ee6d392cd51649c928bd604acafbab4b9d28feb1)]
  + BAM [[Commit Detail](https://github.com/akikuno/DAJIN2/commit/1f3a9812756f0a2607ece3551740e4c67955324c)]

## 📝 Documentation

+ Add a description of the procedure for accepting files generated by Dorado basecaller as input. Issue: #37 [[Commit Detail](https://github.com/akikuno/DAJIN2/commit/c9ebc020fa60980ba7aaaf9295975775ec07da6d)]


## 🔧 Maintenance

+ Specify the Python version to be between 3.8 and 3.10. [[Commit Detail](https://github.com/akikuno/DAJIN2/commit/5fae947eff7da0f7e1ed5e4ff3f95c911fd9f646)]

+ Change `mutation_exporter.report_mutations` to return list[list[str]]. Update the tests accordingly. [[Commit Detail](https://github.com/akikuno/DAJIN2/commit/7153cb143d621e136ca94bfe6b391f1d7b61d438)]

+ Apply formatting with Ruff [[Commit Detail](https://github.com/akikuno/DAJIN2/commit/aec9b697863ef06b4e86e248bebde6616f4eb54e)]

## 🐛 Bug Fixes

+ Add `reallocate_insertion_within_deletion` into `report.mutation_exporter` and reflected it in the mutation info. [[Commit Detail](https://github.com/akikuno/DAJIN2/commit/ed6a96e01bb40c77df9cd3a17a4c29524684b6f1)]

</details>

<details>
<summary> v0.5.0 (2024-06-05) </summary>

## 📝 Documentation

+ Update the issue template from md to yml and modify it to make it easier for users to fill out each item.  [[Commit Detail](https://github.com/akikuno/DAJIN2/commit/ee2c1784e3cb0e72fd09b7c7df577082c19c1a88)]


## 💥 Breaking

+ Extremely low-frequency alleles (less than 0.05%) are considered Nanopore sequence errors and are not clustered #36.
  + Configure `clustering.extract_labels` so that alleles with a low number of reads (0.05% or fewer or 5 reads or fewer) are not clustered. [[Commit Detail](https://github.com/akikuno/DAJIN2/commit/e671b5c84b4cf522faf51823e36fe075b049efcf)]
  + Change `clustering.clustering` to stop if the minimum value of the elements in the cluster is 0.5% or less. [[Commit Detail](https://github.com/akikuno/DAJIN2/commit/74609d5048a4ad8d7004886bf411b2ed4be7fa4b)]
  + Add `consensus.remove_minor_alleles` to remove minor alleles with fewer than 5 reads or less than 0.5% [[Commit Detail](https://github.com/akikuno/DAJIN2/commit/70f675d6a8ea90e9fca51639ddb2b4609e0f4c80)]


+ Save subsetted fastq of a control sample if the read number is too large (> 10,000 reads). The control will have a maximum of 10,000 reads to avoid excessive computational load. [[Commit Detail](https://github.com/akikuno/DAJIN2/commit/d21827f8bbeec326fa2aa4f28feadd6fdecaf554)]

+ If the read length is 500 bases or less, change the mappy preset to `sr`. [[Commit Detail](https://github.com/akikuno/DAJIN2/commit/6e56804ad40780e200f4e9c9ea23294b95443aba)]

+ Update `extract_best_preset` to prioritize `map-ont` and remove `splice` preset if inversion is observed. [[Commit Detail](https://github.com/akikuno/DAJIN2/commit/aa7f6925d6ef4a80a1ba0bbf2b75d8e549ae9863)]


Update the algorithms of `cssplits_hander.reallocate_insertion_within_deletion` to automate change point detection by incorporating temporal changes. [[Commit Detail](https://github.com/akikuno/DAJIN2/commit/7ed8ac8404d18b86c163c71ded6dd1ba784bce79)]



## 🔧 Maintenance

+ Update `deploy_pypi.yml` to use the latest version of Actions. Refer to [the latest official YAML for guidance](https://docs.github.com/actions/automating-builds-and-tests/building-and-testing-python#publishing-to-package-registries). [[Commit Detail](https://github.com/akikuno/DAJIN2/commit/1a54b40146acd21eee30a3a373c44b419d170ad4)]


+ Integrate `requirements.txt` and `MANIFEST.in` into `pyproject.toml` by replacing `setup.py` [[Commit Detail](https://github.com/akikuno/DAJIN2/commit/12f255c3a280098f0310755c51e966031c724932)]

+ Modify to record the execution command of DAJIN2 in the log file [[Commit Detail](https://github.com/akikuno/DAJIN2/commit/38c97a725f6dd3f00162325bf504142f8f8d6594)]

+ Add a test to check if the version in `test_version.sh` matches the version in `pyproject.toml` and `utils.config` [[Commit Detail](https://github.com/akikuno/DAJIN2/commit/a06cb4593ef1a11b3c9826f7ca5532a1bf83f67f)]


+ Rename `consensus.subset_clust` to `consensus.downsample_by_label` to clarify the function's purpose. [[Commit Detail](https://github.com/akikuno/DAJIN2/commit/f6e3f0bc2982996a7dbbc4126a80a7dedd076430)]


+ Update `extract_unique_insertions` to merge highly similar extracted insertion sequences. [[Commit Detail](https://github.com/akikuno/DAJIN2/commit/50fe99f42bcd0bae85bcd0eb4ee371a65f38ea14)]
  + Fix `extract_unique_insertions`: There is a bug where removing the key twice in fasta_insertions_unique caused the index and key to become misaligned in enumerate(distances) if i != key. Therefore, the removal of keys from fasta_insertions_unique is now done all at once at the end. [[Commit Detail](https://github.com/akikuno/DAJIN2/commit/162f248b4deee8c35512b84ec428baec65fd8466)]


+ Add control characters for `fastx_handler.sanitize_filename` as forbidden chars. [[Commit Detail](https://github.com/akikuno/DAJIN2/commit/6b74fce0caa0580c4629a132406206c27a66274d)]


+ Changed the naming convention for the temporary directory: `<sample_name>/<process_content>/<allele_name>/(<label_name>)/file_name`. Example: `flox/consensus/control/1/mutation_loci.pickle`. [[Commit Detail](https://github.com/akikuno/DAJIN2/commit/54fee2f48564c6a29fd5c4151126ba4246e9547c)]

+ Move `sanitze_name` function from `utils.fastx_handler`to `utils.io` [[Commit Detail](https://github.com/akikuno/DAJIN2/commit/a78bd5c0ad8f26bafe369da69607faf9a467c039)]


## 🐛 Bug Fixes

+ Removed `sam_handler.remove_overlapped_reads` to prevent unnecessary trimming of reads. [[Commit Detail](https://github.com/akikuno/DAJIN2/commit/a8991edc0620412c384760d0862e34cc4ea6c0f1)]

+ Fix `preprocess.insertions_to_fasta.remove_minor_groups` to delete the keys (insertion loci) when insertions are removed and result in an empty dict. This prevents errors when accessing non-existent keys in `subset_insertions`. [[Commit Detail](https://github.com/akikuno/DAJIN2/commit/ae8d887282035552c8fbe5c587e43844d5199952)]

+ Fix the bug in `cssplits_handler.convert_cssplits_to_cstag` where the insertion cs tag is not merged with the next cs tag if they have the same operator (e.g., `+A|+A|=T, =T`: before: `+aa=T=T`, after: `+aa=TT`). [[Commit Detail](https://github.com/akikuno/DAJIN2/commit/02d1b4c128004e02671e833136508328d699f53f)]

+ Modified the system to separate intermediate files using a directory structure instead of underscores (`_`), ensuring that no errors occur even if users use allele names containing underscores [[Commit Detail](https://github.com/akikuno/DAJIN2/commit/f70948315a114e0c182895ba4320233f26fc1025)]
  + Thank you @geedrn for reporting the issue #39!

</details>

<details>
<summary> v0.4.6 (2024-05-17) </summary>

## 💥 Breaking

+ Update the log file [Commit Detail](https://github.com/akikuno/DAJIN2/commit/f179c264193391e27f16c66d0f0153f8ae366005)
  + Add the version of DAJIN2 to the log file to track the version of the analysis.
  + Rename the log file to `DAJIN2_log_<current time>.txt` from `<current time>_DAJIN2.log` to enabling open the file in any text editor.

+ Update `mutation_extractor.is_dissimilar_loci` [Commit Detail](https://github.com/akikuno/DAJIN2/commit/2e141bfbbf41a8fe72d11acf159e1974143b7f4e)
  - Rename to `is_dissimilar_loci` from `identify_dissimilar_loci` to explicitly indicate that a boolean is returned.
  - Changed to use cosine distance instead of cosine similarity to make "difference from control" more intuitive.
  - Added a condition to ensure that the cosine distance is not dependent on the specific index: Calculate the cosine distance for 10 bases starting from the neighbor of the corresponding indel, and add the condition that the cosine distances of these adjacent 10 bases should be similar.

+ Update `preprocess.insertions_to_fasta.py` which detects unintended insertion alleles. [Commit Detail](https://github.com/akikuno/DAJIN2/commit/d8bbd9f50b163b6099a5e77c9f7f4de2f5fc08f7)
  + `clustering_insertions`: To accelerate MeanShift clustering, set `bin_seeding=True`. Additionally, because clustering decoys without variation becomes extremely slow, we have switched to using decoys that include slight variations.
  + `extract_unique_insertions`: Within `unintended insertion alleles`, alleles similar to the `intended allele` provided by the user are now excluded.
    + The similarity is defined as there being differences of more than 10 bases

+ Update `preprocess.insertions_to_fasta.clustering_insertions` to consider the length of each insertion sequence during clustering. This allows two alleles, such as `N,(30-base Insertion)` and `(30-base Insertion),N`, to be weighted with different scores as [(1, 30), (30, 1)], enabling correct clustering. [Commit Detail](https://github.com/akikuno/DAJIN2/commit/d41617d8386aa2a4f057cf44c293a1097fa146b6)

+ Update `preprocess.homopolymer_handler`: Scaling data to [0, 1] for cosine similarity, normalizing to match scales due to significant differences in mutation rates between samples and controls. [Commit Detail](https://github.com/akikuno/DAJIN2/commit/0ad27ca2fa7a12ce0cb80e938bc55c903113018f)

## 📝 Documentation

+ Add the descriptions about required Python version supporting from 3.8 to 3.10 due to a Bioconda issue to the README.md. [Commit Detail](https://github.com/akikuno/DAJIN2/commit/0b2f9bf8354e7ff72cc8f8925e1cae6dfba67468)


+ Enhance the descriptions in GitHub Issue templates to clarify their purpose. [Commit Detail](https://github.com/akikuno/DAJIN2/commit/08f3c71bf9f8b755e718eea79dd4a2562aa59297)


## 🔧 Maintenance

+ Move `DAJIN2_VERSION` to `utils.config.py` from `main.py` to make it easier to recognize its location. [Commit Detail](https://github.com/akikuno/DAJIN2/commit/f179c264193391e27f16c66d0f0153f8ae366005)

+ Update `io.read_csv` to return a `list[dict[str, str]]`, not `list[str]` to align the output format with `read_xlsx`. [Commit Detail](https://github.com/akikuno/DAJIN2/commit/d406d34fe990776b6dcecc306ba6fb521c9d0ea0)

+ Update `utils.input_validator` and `preprocess.genome_fetcher` to temporarily disable SSL certificate verification, allowing access to UCSC servers. [Commit Detail](https://github.com/akikuno/DAJIN2/commit/0392fb3fd5c7b87a0773c249ea6e496f69c5af35)

+ Add an example of flox knockin design to the `examples` [Commit Detail](https://github.com/akikuno/DAJIN2/commit/972c3e1b0d9cf04f9ff0d07dd0aaf29deef3b814)


+ Update `preprocess.insertions_to_fasta.py`: The label names for the insertions were not starting from 1, so they have been revised to begin at 1. [Commit Detail](https://github.com/akikuno/DAJIN2/commit/64721e353983447450357b26e0ce5b1ff949d865)

+ Change installer from pip to conda to install mappy in macos-latest (macos-14-arm64) in Github Action [Commit Detail](https://github.com/akikuno/DAJIN2/commit/e1bf83d8f356b5ab5144501de432f83a8394fb16)

## 🚀 Performance

+ Update `consensus.similarity_searcher` to cache onehot encoded controls to avoid redundant computations and increase processing speed. [Commit Detail](https://github.com/akikuno/DAJIN2/commit/0f96c69099cf97e2f4f5a795e224a887c4c667f9)

## 🐛 Bug Fixes

+ Debug `clustering.strand_bias_handler` [Commit Detail](https://github.com/akikuno/DAJIN2/commit/33e955f4afbfa5c30e3494ee97d1fffe33769778)
  + For `positive_strand_counts_by_labels: dict`, there was a bug that caused an error and halted execution when accessing a non-existent key. It has been fixed to output 0 instead.
  + Created a wrapper function `annotate_strand_bias_by_labels` for outputting strand bias. Fixed a bug where the second and subsequent arguments were not being correctly passed when reallocating clusters with strand bias.

+ Fix `preprocess.knockin_handler` to correctly identify the flox knock-in sites as deletions not present in the control.  [Commit Detail](https://github.com/akikuno/DAJIN2/commit/d4d267c99f8c51d3a3f88f67882bead66685f710)

+ Bug fix and update `reallocate_insertion_within_deletion` [Commit Detail](https://github.com/akikuno/DAJIN2/commit/2f356546999f645a8cb8d33a1fc2f64bc6742113)
  - In the script that considers the region between two deletions as an insertion sequence, the size of the other deletion was not taken into account. Even if there was a single base deletion, the entire sequence between the deletions was considered as an insertion sequence. 
  - Therefore, the region between two deletions is now defined as (1) identifying bins where deletions are enriched within appropriate bins (500 bp) continuously, and (2) extracting the precise break points from the start and end of these bins, implementing an algorithm to extract the large deletion region.


</details>

<!--  ------------------------------------------------------------- -->

<details>
<summary> v0.4.5 (2024-04-24) </summary>

## 🐛 Bug Fixes

+ In version 0.4.4 of strand_bias_handler.remove_biased_clusters, there was an error in the continuation condition for removing biased clusters, which has now been corrected. The correct condition should be 'there are alleles with and without strand bias **and** the iteration count is less than or equal to 1000'. Instead, it was incorrectly set to 'there are alleles with and without strand bias **or** the iteration count is less than or equal to 1000'. [Commit Detail](https://github.com/akikuno/DAJIN2/commit/b72b3855121d0da6ac80636089315ecc26464657)

</details>


<details>
<summary> v0.4.4 (2024-04-23) </summary>

## 💥 Breaking

+ Update the threshold from 5 to 0.5 at `identify_dissimilar_loci` to capture 1% minor alleles. [Commit Detail](https://github.com/akikuno/DAJIN2/commit/257b63819921dcf822b831d733f556acd4fec718)

+ Return smaller allele clustering labels (`labels_previous`) when the adjusted Rand index is sufficiently high to reduce predicted allele numbers.
 [Commit Detail](https://github.com/akikuno/DAJIN2/commit/8872daad03bc76acc80fb79fa7260dba73186fae)

## 🔧 Maintenance

+ Add the detailed discription at `identify_dissimilar_loci` to clarify the purpose of the function. [Commit Detail](https://github.com/akikuno/DAJIN2/commit/d2309a133e1bd1f09366477c830923b20e10ca6a)

+ Update a function name of `utils.io.check_excel_or_csv` to `utils.io.determine_file_type` for clarity. [Commit Detail](https://github.com/akikuno/DAJIN2/commit/38f3e2f429eadb3f16dc5f0f64e9b5b135d2cac0)

+ Update examples: In tyr_c230gt_01, the point mutation of Tyr was previously 0.7%, but has been increased to 1.0% by adding point mutation reads from tyr_c230gt_50. [Commit Detail](https://github.com/akikuno/DAJIN2/commit/8967dfd9cc79679be8c7a3e1052467bc57cc375b)

+ Rename `validate_columns_of_batch_file` in test_main.py. [Commit Detail](https://github.com/akikuno/DAJIN2/commit/fc7dc3b9799831b17753f5bbfbd3ca0b4d99e454)

+ Add tests of `strand_bias_handler` [Commit Detail](https://github.com/akikuno/DAJIN2/commit/50703a253f6fde01a002909a3f484141363bbab5)

+ Add type hints and comments in `return_labels` [Commit Detail](https://github.com/akikuno/DAJIN2/commit/02fd72d040865c1c8c81015d965ae12f6788b422)


</details>

<details>
<summary> v0.4.3 (2024-03-29) </summary>

<!-- ## 💥 Breaking -->
## 📝 Documentation

+ Update example dataset and a description of README.md/README_JP.md [Commit Detail](https://github.com/akikuno/DAJIN2/commit/2f9b57057f978b7870e80179c035564c4ee54a40)


<!-- ## 🚀 New Features -->
## 🐛 Bug Fixes

+ Update `preprocess.genome_fetcher_fetch_seq_coordinates` to accurately verify that the entire length of the input sequence is present within the reference sequence. Previously, partial 100% matches were inadvertently accepted; this revision aims to ensure the full alignment of the input sequence with the reference. [Commit Detail](https://github.com/akikuno/DAJIN2/commit/25584734e21e2c8da92d1de12bce498dfc341d03)

+ Update `report.bam_exporter` to be case-sensitive and consistent with directory names. This is to avoid errors caused by the difference between report/bam and report/BAM on Ubuntu, which is case-sensitive to directory names. [Commit Detail](https://github.com/akikuno/DAJIN2/commit/011b21ab32b6965a65e9b442bbf3f2854a44db8e)
  + Thank you @takeiga for reporting the issue #24 !


## 🔧 Maintenance

+ Change `threshold_readnumber` at `labem_merger.merge_labels` from 10 to 5 to capture 1% alleles from 500 total reads. [Commit Detail](https://github.com/akikuno/DAJIN2/commit/8448a8ec1f9efd4d15687a695ab993dc0a27efae)

+ Update the `requirements.txt` to install a newer version of the library. [Commit Detail](https://github.com/akikuno/DAJIN2/commit/d1cbf95b6a16ea720e0033e9a125d6201b99bcee)

+ Update `report.report_bam` and rename to `report.bam_exporter`: [Commit Detail](https://github.com/akikuno/DAJIN2/commit/011b21ab32b6965a65e9b442bbf3f2854a44db8e)
  + Use UUID instead of random number for the temporary file name.
  + Rename `realign` to `recalculate_sam_coodinates_to_reference` for the readability of the function name.
  + Add `convert_pos_to_one_indexed` to convert the 0-based position to 1-based position and suppress samtools warning.
    + Warning: `[W::sam_parse1] mapped query cannot have zero coordinate; treated as unmapped`
  + Add tests for the `write_sam_to_bam` function

+ Move `read_sam` function from sam_handler to io module. [Commit Detail](https://github.com/akikuno/DAJIN2/commit/f9b9382ab706530b0cd4c34d7ff8f8c79002b654)

+ Rename `report.report_mutation`, `report.report_files` to `report.mutation_exporter` and `report.sequence_exporter` to be more explicit. [Commit Detail](https://github.com/akikuno/DAJIN2/commit/35d8250876cd845623e63c898d7c608d27a82a45)

</details>


<details>
<summary> v0.4.2 (2024-03-25) </summary>

## 🔧 Maintenance

+ Remove multi-mapping reads, as multi-mapping reads are mostly reads that are locally mapped to low-complexity regions. [Commit Detail](https://github.com/akikuno/DAJIN2/commit/d00bea83366113ff0ccf80639b75bb7edbb4ed2f)

+ Create `preprocess.input_formatter.py` to summarize formatting functions to a module. [Commit Detail](https://github.com/akikuno/DAJIN2/commit/bb45bb81a8deb530109de18e794f63ecb088f651)

+ Refactor `directory_manager.py` [Commit Detail](https://github.com/akikuno/DAJIN2/commit/9d558304609935f9d3320cc1f6d7b3a46168d9e2)

+ Refactor `preprocess.__init__.py` [Commit Detail](https://github.com/akikuno/DAJIN2/commit/43ab68a135498b3e8192e1facbd085152e429f86)

+ To increase cohesion by functions of the same category into a single module, we have migrated `preprocess.fastx_parser` to `utils.fastx_handler`. [Commit Detail](https://github.com/akikuno/DAJIN2/commit/e9396369c47cb09af7d78c0f9eb71a5f225232e5)

+ Remove the packages that are no longer in use from `requirements.txt`. [Commit Detail](https://github.com/akikuno/DAJIN2/commit/43ab68a135498b3e8192e1facbd085152e429f86)

+ Add `read_sam` in sam_handler module. [Commit Detail](https://github.com/akikuno/DAJIN2/commit/b37b3750f76ef354827229a7467e56a439225fe1)

+ Revise the docstring of `export_fasta_files`. [Commit Detail](https://github.com/akikuno/DAJIN2/commit/4c6fa03f61d8473e50c187e3bf4cb3e8685f2631)

+ Standardize to use `dataclass` instead of `NamedTuple`. [Commit Detail](https://github.com/akikuno/DAJIN2/commit/b7c34fbcda51ef037488f1f58564fa72128033f1)

</details>


<details>
<summary> v0.4.1 (2024-02-13) </summary>

## 📝 Documentation

- Added documentation for a new feature in `README.md`: DAJIN2 can now detect complex mutations characteristic of genome editing, such as insertions occurring in regions where deletions have occurred.

## 🚀 New Features

- Introduced `cssplits_handler.detect_insertion_within_deletion` to extract insertion sequences within deletions. This addresses cases where minimap2 may align bases that partially match the reference through local alignment, potentially failing to detect them as insertions. This enhancement ensures the proper detection of insertion sequences. [Commit Detail](https://github.com/akikuno/DAJIN2/commit/7651e20852b94ed4d5bb38539bb56229dcc8b763)

- Added `report.insertion_refractor.py` to include original insertion information in the consensus for mappings made by insertion. This addition enables the listing of both insertions and deletions within the insertion allele on a single HTML file. [Commit Detail](https://github.com/akikuno/DAJIN2/commit/e6c3b636bb2ba537d1341d1042341afd6583dd0b)

## 🔧 Maintenance

- Updated `insertions_to_fasta.py`. [Commit Detail](https://github.com/akikuno/DAJIN2/commit/7927feb0bb4f3091537aaebabd60a441456a3413)
  - Modified the approach to reduce randomness by replacing set or frozenset with list or tuple, and using `random.sample()` for subsetting reads.
  - Refactored `call_consensus_insertion_sequence`.
  - Fixed a bug in `extract_score_and_sequence` to ensure correct appending of scores for the insertions_merged_subset.

- Changed the function name of `report` to be more explicit. [Commit Detail](https://github.com/akikuno/DAJIN2/commit/93132c5beba17278c7d67b76817bb13dfaae57a3)

- Updated `utils.report_report_generator` [Commit Detail](https://github.com/akikuno/DAJIN2/commit/821f06f05b5ed2f4ba2d7baad6159d774d2e5db0)
  - Capitalized "Allele" (e.g., control) and "Allele type" (e.g., intact).
  - Changed the output format of read_all and read_summary from CSV to XLSX.
  - Corrected the order of the Legend to follow a logical sequence from control to sample, and then to specific insertions.

- Updated `utils.io.read_xlsx` to switch from using pandas to openpyxl due to the DeprecationWarning in Pandas being cumbersome. [Commit Detail](https://github.com/akikuno/DAJIN2/commit/5d942bace8417bb973441b360a0ec31d77d81e24)

## 🐛 Bug Fixes

- Added `=` to the prefix for valid cstag recognition when there is an `n` in inversion. [Commit Detail](https://github.com/akikuno/DAJIN2/commit/747ff3ece221a8c1e4f1ba1b696c4751618b4992)

- Modified the io.load_from_csv function to trim spaces before and after each field, addressing an error caused by spaces in batch.csv. [Commit Detail](https://github.com/akikuno/DAJIN2/commit/f5d49230f8ebd37061a27d6767d3c1954b8f8576)

## ⛔️ Deprecated

- Removed `reads_all.csv`. This CSV file, which showed the allele for each read, is no longer reported due to its limited usefulness and because the same information can be obtained from the BAM file. [Commit Detail](https://github.com/akikuno/DAJIN2/commit/76e3eaee320deb79cbf3cf97cc6aed69c5bbc3ef)

</details>


<details>
<summary> v0.4.0 (2024-01-20) </summary>

## 💥 Breaking

+ Changed the input from a path to a FASTQ file to **a path to a directory**: The output of Guppy is now stored in multiple FASTQ files under the `barcodeXX/` directory. Previously, it was necessary to combine the FASTQ files in the `barcodeXX/` directory into one and specify it as an argument. With this revision, it is now possible to directly specify the `barcodeXX` directory, allowing users to seamlessly proceed to DAJIN2 analysis after Guppy processing.
[Commit Detail](https://github.com/akikuno/DAJIN2/commit/d35ce6f89278d0361cc2b5b30fecfabbc66aa1c4)

## 📝 Documentation

+ Changed `conda config --set channel_priority strict` to `conda config --set channel_priority flexible` for installation process in TROUBLESHOOTING.md. [Commit Detail](https://github.com/akikuno/DAJIN2/commit/c95681a8f2b6e725b0b737498981ad767eab842c)

## 🚀 New Features

+ Apple Silicon (ARM64) supoorts. [Commit Detail](https://github.com/akikuno/DAJIN2/commit/435bab6c56cb2172601d4b37488850fe48046f9c)

+ Changed the definition of the minor allele from a read number of less than or equal to 10 to less than or equal to 5. This is based on the assumption that one sample contains 1000 reads, where 0.5% corresponds to 5 reads. [Commit Detail](https://github.com/akikuno/DAJIN2/commit/80a3ddcf7cac3eed2bcc76b88ea534873af4dd90)


## 🔧 Update

+ Update `preprocess.insertion_to_fasta` to facilitate the discrimination of Insertion alleles, the Reference for Insertion alleles has been saved in FASTA/HTML directory. [Commit Detail](https://github.com/akikuno/DAJIN2/commit/5899543077f0398863b6316d8c3e953b5f125f55)

+ Update `insertions_to_fasta.extract_enriched_insertions`: Previously, it calculated the presence ratio of insertion alleles separately for samples and controls, filtering at 0.5%. However, due to a threshold issue, some control insertions were narrowly missing the threshold, resulting in them being incorrectly identified as sample-specific insertions. To rectify this, the algorithm now clusters samples and controls together, excluding clusters where both types are mixed. This modification allows for the extraction of sample-specific insertion alleles. [Commit Detail](https://github.com/akikuno/DAJIN2/commit/65030daba7c56a6c3f3f685832084b71c6b2e1c3)

+ Updated `preprocess.insertions_to_fasta.count_insertions` of the counting method to treat similar insertions as identical. Previously, the same insertion was erroneously counted as different ones due to sequence errors. [Commit Detail](https://github.com/akikuno/DAJIN2/commit/7bc18f486253e876d51a296f64909e1c73114e79)

+ Updated `preprocess.insertions_to_fasta.merge_similar_insertions`: Previously, clustering was done using MiniBatchKMeans, but this method had an issue where it excessively clustered when only highly similar insertion sequences existed. Therefore, a strategy similar to `extract_enriched_insertions` was adopted, changing the algorithm to one that mixes with a uniform distribution of random scores before clustering. [Commit Detail](https://github.com/akikuno/DAJIN2/commit/fb7074cab9d9e4e3d293cb5487a3525a5faf06fd)

+ Added `preprocess.insertions_to_fasta.clustering_insertions`: Combined the clustering methods used in `extract_enriched_insertions` and `merge_similar_insertions` into a common function. [Commit Detail](https://github.com/akikuno/DAJIN2/commit/6d7ff79351c5f60320b2269accb0e3bc159fdd5b)


+ Moved the `call_sequence` function to the `cssplits_handler` module. [Commit Detail](https://github.com/akikuno/DAJIN2/commit/ef5b0bf41ab33a7e8d06d33fe7fa6c27a443742a)

## 🐛 Bug Fixes

+ Debug `clustering.merge_labels` to be able to correctly revert minor labels back to parent labels. [Commit Detail](https://github.com/akikuno/DAJIN2/commit/8127a94e042328b87e456d3748ebea66a845ba1a)


+ Updated `utils.input_validator.validate_genome_and_fetch_urls` to obtain `available_server` more explicitly. Previously, it relied on HTTP response codes, but there were instances where the UCSC Genome Browser showed a normal (200) response while internally being in error. Therefore, with this change, a more explicit method is employed by searching for specific keywords present in the normal HTML, to determine if the server is functioning correctly. [Commit Detail](https://github.com/akikuno/DAJIN2/commit/24a02591e8a146030012dbf564e4b6cd98d42139)

+ Added `config.reset_logging` to reset the logging configuration. Previously, when batch processing multiple experiment IDs (names), a bug existed where the log settings from previous experiments remained, and the log file name was not updated. However, with this change, log files are now created for each experiment ID. [Commit Detail](https://github.com/akikuno/DAJIN2/commit/b83669c627710a5e358f934212e961373203ee52)

+ Debugged `core.py`: Modified the specification of `paths_predefined_fasta` to accept input from user-entered ALLELE data. Previously, it accepted fasta files stored in the fasta directory. However, this approach had a bug where fasta files left over from a previously aborted run (which included newly created insertions) were treated as predefined. This resulted in new insertions being incorrectly categorized as predefined. [Commit Detail](https://github.com/akikuno/DAJIN2/commit/6dd9247f010eb6168157ae9236a634efcfb84a5f)

</details>

<details>
<summary> v0.3.6 (2024-01-10) </summary>

## 📝 Documentation

- Added a quick guide for installation to TROUBLESHOOTING.md. [Commit Detail](https://github.com/akikuno/DAJIN2/commit/cefed0ff4d04282b9915486be07de85b2b77b657)

## 🚀 Update

### Preprocess

- Updated `input_validator.py`: The UCSC Blat server sometimes returns a 200 HTTP status code even when an error occurs. In such cases, "Very Early Error" is indicated in the title. Therefore, we have made it so that it returns False in those situations. [Commit Detail](https://github.com/akikuno/DAJIN2/commit/4ad9c9ef8bd963a6e20c1721480aed0fe7922760)

- Simplified `homopolymer_handler.py` for error detection using cosine similarity. [Commit Detail](https://github.com/akikuno/DAJIN2/commit/21c2596805c36074f360285600e60ee76b948908)

- Updated `mutation_extractor.py` to use cosine similarity to filter dissimilar loci. [Commit Detail](https://github.com/akikuno/DAJIN2/commit/c9f5aa7b48581e58d99fe8c31275c422756aa9f1)

- Updated the `mutation_extractor.identify_dissimilar_loci` so that it unconditionally returns True if the 'sample' shows more than 5% variation compared to the 'control'. [Commit Detail](https://github.com/akikuno/DAJIN2/commit/0cbec5217fdfba6886979eb86cf970b587e83e5f)

- Added `preprocess.midsv_caller.convert_consecutive_indels_to_match`: Due to alignment errors, instances where a true match is mistakenly replaced with "insertion following a deletion" are corrected. For example, "=C,=T" mistakenly replaced by "-C,+C|=T" is reverted back to "=C,=T". [Commit Detail](https://github.com/akikuno/DAJIN2/commit/69c56fa904ef847dc5b0e2dcdb90303409412d0f)

### Classification

- Added `allele_merger.merge_minor_alleles` to reclassify alleles with fewer than 10 reads to suppress excessive subdivision of alleles. [Commit Detail](https://github.com/akikuno/DAJIN2/commit/b0752960def313e237ccf7d44542f9810cad0c00)

### Clustering

- Added the function `merge_minor_cluster` to revert labels clustered with fewer than 10 reads back to the previous labels to suppress excessive subdivision of alleles. [Commit Detail](https://github.com/akikuno/DAJIN2/commit/4bd9f7dd806d192475d8d4f20c1e50c37281d64e)

- Updated `generate_mutation_kmers` to consider indices not registered in mutation_loci as mutations by replacing them with "@". For example, "=G,=C,-C" and "=G,=G,=C" become "@,@,@" in both cases, making them the same and ensuring they do not affect clustering. [Commit Detail](https://github.com/akikuno/DAJIN2/commit/9eefaaa1a9be3922b60655292c0a310e0f5fc76d)

### Consensus

- Implemented `LocalOutlierFactor` to filter abnormal control reads. [Commit Detail](https://github.com/akikuno/DAJIN2/commit/94c48da01fc039902c97a23ddea47dd5f2b42ab4)

</details>


<details>
<summary> v0.3.5 (2023-12-23) </summary>

## 📝 Documentation

+ [x] Added `ROADMAP.md` to track the progress of the project [Commit Detail](https://github.com/akikuno/DAJIN2/commit/cf05d3e5c9b1d3ee806d66c9c1d9f8079863e312)
+ [x] Added *Prerequisites* section to README.md [Commit Detail](https://github.com/akikuno/DAJIN2/commit/7d5a3cd8305f9d414a492f5223d5dbec7399aa46)

## 🚀 Features

### Preprocessing

+ [x] Updated `homopolymer_handler.get_counts_homopolymer` to change to count mutations in homopolymer regions considering only the control [Commit Detail](https://github.com/akikuno/DAJIN2/commit/e5d061750c66bdc225fcddfae6e2d2a12fe49ad2)

### Clustering

+ [x] Changed clustering algorithm from KMeans to BisectingKMeans to handle larger dataset [Commit Detail](https://github.com/akikuno/DAJIN2/commit/7733524625de77c814496791a461eb7bbff54d0e)

### Consensus

+ [x] Added `convert_consecutive_indels_to_match` to offset the effect when the same base insertion/deletion occurs consecutively [Commit Detail](https://github.com/akikuno/DAJIN2/commit/a678615b4ffeeefdc9509f49651698281b1aff22)

+ [x] Added `similarity_searcher.py` to extract control reads resembling the consensus sequence, thereby enhancing the accuracy of detecting sample-specific mutations. [Commit Detail](https://github.com/akikuno/DAJIN2/commit/98a8a45e13835502f7dea2622274da81bbbc3ba3)

+ [x] Changed the method in `clust_formatter.get_thresholds`` to dynamically define the thresholds for ignoring mutations, instead of using fixed values.[Commit Detail](https://github.com/akikuno/DAJIN2/commit/2249d1601ad619a7db0fcc9ebf79d63f8dcf164b)

+ [x] Removed code that was previously commented out [Commit Detail](https://github.com/akikuno/DAJIN2/commit/2249d1601ad619a7db0fcc9ebf79d63f8dcf164b)

+ [x] Add `is_consensus` argument: When it comes to consensus, if the difference between sample and control is more than 20%, it is unconditionally considered a mutation. [Commit Detail](https://github.com/akikuno/DAJIN2/commit/7bca4590f97e1858304c3e9fb66c54a279dfcdf0)


## 🐛 Bug Fixes

+ None

## 🔧 Maintenance

+ [x] Modified batch processing to run on a single CPU thread per process [Commit Detail](https://github.com/akikuno/DAJIN2/commit/7b43e36b9482cceabe79f47814f62f69d46b7d3e)

+ [x] Simplifed import path [Commit Detail](https://github.com/akikuno/DAJIN2/commit/6e2d1726edc49fc638b87526a3f4fcbf1eead4e0)
  + `preprocess.midsv_caller.execute` to `preprocess.generate_midsv`
  + `preprocess.mapping.generate_sam` to `preprocess.generate_sam`

+ [x] Added tests to `consensus.convert_consecutive_indels_to_match` [Commit Detail](https://github.com/akikuno/DAJIN2/commit/c4932dc1c0776b604122558331a9fb41a29244af)

## ⛔️ Deprecated

+ None

</details>
