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
- barcode09 allele1 の`N`
- barcode11 allele2 の大型欠失が反映されていない
- barcode28 allele1 の`N`
- FASTQ、VCFを出力する
 -->
# v0.3.6 (yyyy-mm-dd)

## 📝 Documentation

+ Added a quick quide of installation to TROUBLESHOOTING.md [Commit Detail](https://github.com/akikuno/DAJIN2/commit/cefed0ff4d04282b9915486be07de85b2b77b657)

## 🚀 Features

### Preprocess

+ Simplyfy `homopolymer_handler.py` for the error detection using cosine similarity [Commit Detail](https://github.com/akikuno/DAJIN2/commit/21c2596805c36074f360285600e60ee76b948908)

+ Use cosine similarity to filter dissimilar loci at `mutation_extractor.py` [Commit Detail](https://github.com/akikuno/DAJIN2/commit/c9f5aa7b48581e58d99fe8c31275c422756aa9f1)

### Classification

#### Skipped the classification of minor alleles to suppress excessive subdivision of alleles.

+ Added `merge_minor_alleles` to reclassify alleles with less than 10 reads. [Commit Detail](https://github.com/akikuno/DAJIN2/commit/a9399f5d4aaa0f6bcdbea0ead39b2aabf8223f5d)

### Clustering

#### Skipped the clustering of minor alleles to suppress excessive subdivision of alleles.

+ Added the function `merge_minor_cluster` to revert labels clustered with less than 10 reads back to the previous labels.
 [Commit Detail](https://github.com/akikuno/DAJIN2/commit/4bd9f7dd806d192475d8d4f20c1e50c37281d64e)

### Consensus

+ Use `LocalOutlierFactor` to filter abnormal control reads [Commit Detail](https://github.com/akikuno/DAJIN2/commit/4bd9f7dd806d192494c48da01fc039902c97a23ddea47dd5f2b42ab475d8d4f20c1e50c37281d64e)


## 🐛 Bug Fixes

### Consensus

+ 大型欠失の内部で欠失が反映されないバグを修正 [Commit Detail](https://github.com/akikuno/DAJIN2/commit/XXX)

## 🔧 Maintenance


## ⛔️ Deprecated

---

# 💡 Future Tasks

+ Remove minor alleles with predicted insertion
+ Enhance the Clarity of Insertion Allele Identification.
+ Develop and Integrate Inversion Detection Capability

-------------

# Past Logs

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
