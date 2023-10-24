[![License](https://img.shields.io/badge/License-MIT-9cf.svg?style=flat-square)](https://choosealicense.com/licenses/mit/)
[![Test](https://img.shields.io/github/actions/workflow/status/akikuno/dajin2/pytest.yml?branch=main&label=Test&color=brightgreen&style=flat-square)](https://github.com/akikuno/dajin2/actions)
[![Python](https://img.shields.io/pypi/pyversions/DAJIN2.svg?label=Python&color=blue&style=flat-square)](https://pypi.org/project/DAJIN2/)
[![PyPI](https://img.shields.io/pypi/v/DAJIN2.svg?label=PyPI&color=orange&style=flat-square)](https://pypi.org/project/DAJIN2/)
[![Bioconda](https://img.shields.io/conda/v/bioconda/dajin2?label=Bioconda&color=orange&style=flat-square)](https://anaconda.org/bioconda/dajin2)
[![DOI](https://zenodo.org/badge/387721337.svg)](https://zenodo.org/badge/latestdoi/387721337)


<p align="center">
<img src="https://user-images.githubusercontent.com/15861316/261833016-7f356960-88cf-4574-87e2-36162b174340.png" width="90%">
</p>

DAJIN2は、ゲノム編集技術で作製された動物の遺伝型を解析するためのソフトウェアです。

## 🌟 特徴

- ゲノム編集イベントを網羅的に検出することができます。  
- ゲノム編集結果を可視化し、直観的に確認することができます。  
- 多サンプル処理が可能です。  


## ⚠️ 諸注意

DAJIN2は開発段階のツールです。  
点変異、欠失、挿入(flox)に関する基本的なテストは完了していますが、複雑な編集（例：逆位）には対応していません。  
ご利用中にエラーに遭遇した場合、お手数ですが[Issues](https://github.com/akikuno/DAJIN2/issues)からご報告をお願いします（日本語でも構いません）。  

## 🛠 インストール

### [Bioconda](https://anaconda.org/bioconda/DAJIN2) （推奨）

```bash
conda install -c bioconda DAJIN2
```

### [PyPI](https://pypi.org/project/DAJIN2/)

```bash
pip install DAJIN2
```

> **Warning**  
> pip経由でインストールすると、`Failed to build mappy`というエラーが出ることがあります。  
> 下記のコードで`gcc`と`zlib`をインストールしたのちに、再度DAJIN2のインストールを行ってください。  
> `sudo apt install gcc zlib1g zlib1g-dev` (Ubuntu)  
> `brew install gcc zlib` (macOS)
```


## 💡 使用法

### 単一サンプルの解析

単一サンプル（サンプルのFASTQとコントロールのFASTQ）の解析手順は以下の通りです。


```bash
DAJIN2 <-s|--sample> <-c|--control> <-a|--allele> <-n|--name> [-g|--genome] [-t|--threads] [-h|--help] [-v|--version]

options:
  -s, --sample              Path to a sample FASTQ file
  -c, --control             Path to a control FASTQ file
  -a, --allele              Path to a FASTA file
  -n, --name                Output directory name
  -g, --genome (Optional)   Reference genome ID (e.g hg38, mm39) [default: '']
  -t, --threads (Optional)  Number of threads [default: 1]
  -h, --help                show this help message
  -v, --version             show the version number
```

#### 使用例

```bash
# Donwload the example dataset
wget https://github.com/akikuno/DAJIN2/raw/main/examples/example-single.tar.gz
tar -xf example-single.tar.gz

# Run DAJIN2
DAJIN2 \
    --name stx2-deletion \
    --sample example-single/sample.fq.gz \
    --control example-single/control.fq.gz \
    --allele example-single/design.fa \
    --genome mm39 \
    --threads 10

# 2023-06-04 11:30:03: example-single/control.fq.gz is now processing...
# 2023-06-04 11:30:06: Preprocess example-single/control.fq.gz...
# 2023-06-04 11:30:06: Mapping example-single/control.fq.gz...
# 2023-06-04 11:30:21: Call MIDSV example-single/control.fq.gz...
# 2023-06-04 11:30:31: 🍵 example-single/control.fq.gz is finished!
# 2023-06-04 11:30:31: example-single/sample.fq.gz is now processing...
# 2023-06-04 11:30:35: Preprocess example-single/sample.fq.gz...
# 2023-06-04 11:34:13: Classify example-single/sample.fq.gz...
# 2023-06-04 11:34:18: Clustering example-single/sample.fq.gz...
# 2023-06-04 11:35:01: Consensus calling example-single/sample.fq.gz...
# 2023-06-04 11:35:08: 🍵 example-single/sample.fq.gz is finished!
# 🎉 Finished! Open DAJIN_Results/stx2-deletion to see the report.
```

### 複数サンプルの一括処理

`batch`サブコマンドを利用することで、複数のFASTQファイルを同時に処理することができます。  
この際、サンプル情報をまとめたCSVファイルやExcelファイルが必要となります。  
具体例としては、[こちら](https://github.com/akikuno/DAJIN2/blob/main/examples/example-batch/batch.csv)を参照してください。


```bash
DAJIN2 batch <-f|--file> [-t|--threads] [-h]

options:
  -f, --file                Path to a CSV or Excel file
  -t, --threads (Optional)  Number of threads [default: 1]
  -h, --help                Show this help message
```

#### 使用例

```bash
# Donwload the example dataset
wget https://github.com/akikuno/DAJIN2/raw/main/examples/example-batch.tar.gz
tar -xf example-batch.tar.gz

# Run DAJIN2
DAJIN2 batch --file example-batch/batch.csv --threads 3

# 2023-07-31 17:01:10: example-batch/tyr_control.fq.gz is now processing...
# 2023-07-31 17:01:16: Preprocess example-batch/tyr_control.fq.gz...
# 2023-07-31 17:01:48: Output BAM files of example-batch/tyr_control.fq.gz...
# 2023-07-31 17:01:52: 🍵 example-batch/tyr_control.fq.gz is finished!
# 2023-07-31 17:01:52: example-batch/tyr_c230gt_50%.fq.gz is now processing...
# 2023-07-31 17:01:52: example-batch/tyr_c230gt_10%.fq.gz is now processing...
# 2023-07-31 17:01:52: example-batch/tyr_c230gt_01%.fq.gz is now processing...
# 2023-07-31 17:01:55: Preprocess example-batch/tyr_c230gt_01%.fq.gz...
# 2023-07-31 17:01:55: Preprocess example-batch/tyr_c230gt_50%.fq.gz...
# 2023-07-31 17:01:55: Preprocess example-batch/tyr_c230gt_10%.fq.gz...
# 2023-07-31 17:02:17: Classify example-batch/tyr_c230gt_50%.fq.gz...
# 2023-07-31 17:02:19: Clustering example-batch/tyr_c230gt_50%.fq.gz...
# 2023-07-31 17:02:34: Classify example-batch/tyr_c230gt_01%.fq.gz...
# 2023-07-31 17:02:35: Classify example-batch/tyr_c230gt_10%.fq.gz...
# 2023-07-31 17:02:39: Clustering example-batch/tyr_c230gt_01%.fq.gz...
# 2023-07-31 17:02:39: Clustering example-batch/tyr_c230gt_10%.fq.gz...
# 2023-07-31 17:02:53: Consensus calling of example-batch/tyr_c230gt_50%.fq.gz...
# 2023-07-31 17:02:59: Output reports of example-batch/tyr_c230gt_50%.fq.gz...
# 2023-07-31 17:03:04: 🍵 example-batch/tyr_c230gt_50%.fq.gz is finished!
# 2023-07-31 17:03:39: Consensus calling of example-batch/tyr_c230gt_01%.fq.gz...
# 2023-07-31 17:03:51: Output reports of example-batch/tyr_c230gt_01%.fq.gz...
# 2023-07-31 17:04:03: 🍵 example-batch/tyr_c230gt_01%.fq.gz is finished!
# 2023-07-31 17:04:08: Consensus calling of example-batch/tyr_c230gt_10%.fq.gz...
# 2023-07-31 17:04:16: Output reports of example-batch/tyr_c230gt_10%.fq.gz...
# 2023-07-31 17:04:24: 🍵 example-batch/tyr_c230gt_10%.fq.gz is finished!
# 🎉 Finished! Open DAJIN_Results/tyr-substitution to see the report.
```

## 📈 レポート内容

DAJIN2の処理が完了すると、**DAJIN_Results**というディレクトリが作られます。  
このDAJIN_Resultsディレクトリには、以下のファイルが含まれています：  

```text
DAJIN_Results/tyr-substitution
├── BAM
│   ├── tyr_c230gt_01%
│   ├── tyr_c230gt_10%
│   ├── tyr_c230gt_50%
│   └── tyr_control
├── FASTA
│   ├── tyr_c230gt_01%
│   ├── tyr_c230gt_10%
│   └── tyr_c230gt_50%
├── HTML
│   ├── tyr_c230gt_01%
│   ├── tyr_c230gt_10%
│   └── tyr_c230gt_50%
├── MUTATION_INFO
│   ├── tyr_c230gt_01%.csv
│   ├── tyr_c230gt_10%.csv
│   └── tyr_c230gt_50%.csv
├── read_all.csv
├── read_plot.html
├── read_plot.pdf
└── read_summary.csv
```

### 1. BAM

BAMディレクトリには、入力のFASTQと、アレルごとに分類されたreadsのBAMファイルが格納されています。  

> **Note**  
> `genome`オプションで参照ゲノムを指定すると、その参照ゲノムにアライメントされます。  
> 指定がない場合、入力のFASTAファイルのcontrolアレルにアライメントされます。


### 2. FASTA / HTML

FASTAディレクトリには、各アレルのFASTAファイルが保存されます。  
HTMLディレクトリには、変異箇所が色付けされた各アレルのHTMLファイルが保存されます。  
Tyr点変異の例を以下に示します：
- 点変異箇所は**緑色**でハイライトされています。


<img src="https://user-images.githubusercontent.com/15861316/274518501-2ca3f442-1b86-4635-be3d-fd37575c4ca2.png" width="75%" />


### 3. MUTATION_INFO

MUTATION_INFOディレクトリには、各アレルの変異箇所を示すテーブルが保存されます。  
Tyr点変異の例を以下に示します：
- 点変異の染色体上の位置と、変異の種類が記載されています。

<img src="https://user-images.githubusercontent.com/15861316/274519342-a613490d-5dbb-4a27-a2cf-bca0686b30f0.png" width="75%">

### 4. read_plot.html / read_plot.pdf

read_plot.html および read_plot.pdf は、各アレルの割合を図示しています。  
図中の**Allele type**はアレルの種類を、**% of reads**は該当するリードのアレル割合を示しています。  

また、**Allele type**の種類は以下の通りです：

- **intact**：入力のFASTAアレルと完全に一致するアレル
- **indels**：50塩基以内の置換、欠失、挿入、逆位
- **sv**：50塩基以上の置換、欠失、挿入、逆位


<img src="https://user-images.githubusercontent.com/15861316/274521067-4d217251-4c62-4dc9-9c05-7f5377dd3025.png" width="75%">

> **Warning**  
> PCRアンプリコンのシーケンシングでは、増幅バイアスのため、**% of reads**が実際のアレルの割合と一致しないことがあります。  
> 特に大型の欠失が存在する場合、欠失アレルが顕著に増幅されることから、実際のアレル割合を反映していない可能性が高まります。

### 5. read_all.csv / read_summary.csv

- read_all.csv：各リードがどのアレルに分類されたかが記録されています。
- read_summary.csv：各アレルのリード数と存在割合が記述されています。


## 📄 参考文献

[Kuno A, et al. (2022) DAJIN enables multiplex genotyping to simultaneously validate intended and unintended target genome editing outcomes. *PLoS Biology* 20(1): e3001507.](https://doi.org/10.1371/journal.pbio.3001507)
