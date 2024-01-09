# Development Logs of DAJIN2

<!-- TEMPLATE
# v0.0.0 (yyyy-mm-dd)
## 📝 Documentation
## 🚀 Features
## 🐛 Bug Fixes
## 🔧 Maintenance
## ⛔️ Deprecated
+ [ ] XXX [Commit Detail](https://github.com/akikuno/DAJIN2/commit/xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx)
-->

<!-- memo ToDo
- ✅ barcode09 allele1 の`N`
- barcode11 allele2 の大型欠失が反映されていない ＆ barcode28 allele1 の`N` → Control側の問題
- FASTQ、VCFを出力する

+ # 💡 Future Tasks
+ Remove minor alleles with predicted insertion
+ Enhance the Clarity of Insertion Allele Identification.
+ Develop and Integrate Inversion Detection Capability
+ ReferenceのアレルをFASTA/HTMLディレクトリに保存する
 -->




-------------

# Past Logs

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
