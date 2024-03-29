<!-- TEMPLATE
# v0.0.0 (yyyy-mm-dd)
## 💥 Breaking
## 📝 Documentation
## 🚀 New Features
## 🐛 Bug Fixes
## 🔧 Maintenance
## ⛔️ Deprecated
- XXX [Commit Detail](https://github.com/akikuno/DAJIN2/commit/xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx)
-->

<!-- 💡 ToDo
- VCF、PDFを出力する
- 逆位アレルでの検証を加える
- nCATSがほしい…
 -->

# v0.4.2 (2024-03-25)

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


-------------------------------------------------------------

# Past Releases


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
