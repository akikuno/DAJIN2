
<p align="center">
<img src="https://github.com/akikuno/DAJIN/blob/master/misc/images/DAJIN-logo.png" width="90%">
</p>

[![MIT License](http://img.shields.io/badge/license-MIT-blue.svg?style=flat)](LICENSE)

## 特徴

- **移植性**：Windows10 (WSL), Linux, macOS
- **低依存性**：minimap2とsamtoolsのみ必要とします
- **高速**： DAJINの100倍高速です

## セットアップ

### 動作環境

LinuxまたはWindows 10 ([WSL](https://docs.microsoft.com/ja-jp/windows/wsl/install-win10))で動作確認をしています.  
検証済みの環境は[こちら](https://github.com/akikuno/DAJIN/blob/master/docs/TESTED_SYSTEMS.md)です.  

### `conda`をインストールします


wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
chmod +x Miniconda3-latest-Linux-x86_64.sh
./Miniconda3-latest-Linux-x86_64.sh

condaおよびインストールについての詳細は[こちら](https://docs.conda.io/projects/conda/en/latest/)にございます.  


### DAJINをダウンロードします

git clone https://github.com/akikuno/DAJIN.git

または下記URLでZIPファイルをダウンロードしてください.  

https://github.com/akikuno/DAJIN/archive/master.zip

### DAJINをダウンロードします


## 利用方法

### 入力ファイルの用意

以下のような入力ファイルを作製します.  

design=DAJIN/example/design.txt
input_dir=DAJIN/example/demultiplex
control=barcode01
genome=mm10
grna=CCTGTCCAGAGTGGGAGATAGCC,CCACTGCTAGCTGTGGGTAACCC
output_dir=DAJIN_cables2
threads=10
filter=on

各項目の情報は以下のとおりです. 各項目は順不同です.  

- **desing**: 考えられる遺伝型の配列を記載したFASTA形式のテキストファイルです.  ">wt"と">target"の2つは含まれている必要があります. 
- **input_dir**: demultiplex済みのFASTA/FASTQファイルを含むディレクトリです. 
- **control**: 野生型コントロールのバーコード番号です. 
- **genome**: `mm10`, `hg38`等の参照ゲノムです. 
- **grna**: PAMを含むgRNA配列です. 2つ以上の配列はコンマ（,）で区切ります. 
- **output_dir（オプショナル）**: 結果を保存するディレクトリの名前です. デフォルトは`DAJIN_results`です. 
- **threads（オプショナル）**: DAJINに使用するCPUスレッド数です. デフォルトでは`3分の2`を使用します. 
- **filter（オプショナル**: on/off）: マイナーアレル（3%以下）を解析から除きます. デフォルトは"on"です. 


### DAJINの実行

./DAJIN/DAJIN.sh -f [入力ファイルのPATH]

下記のコマンドで例を実行します.

./DAJIN/DAJIN.sh -i DAJIN/example/design.txt

### 結果のレポートについて

DAJINは2つのファイル（`Details.csv`, `Details.pdf`）と2つのフォルダ（`BAM`, `Consensus`）を出力します. 

#### Details.csv

`Details.csv` はアレル情報を記載しています.

| Sample    |  Allele ID |  % of reads |  Allele type  |  Indel |  Large indel |  Design |
|-----------|------------|-------------|---------------|--------|--------------|---------|
| barcode01 | 1          | 100         | wt            | -      | -            | -       |
| barcode02 | 1          | 11.8        | abnormal      | +      | +            | -       |
| barcode02 | 2          | 88.2        | target        | -      | -            | +       |
| barcode03 | 1          | 9.9         | abnormal      | +      | +            | -       |
| barcode03 | 2          | 38.5        | abnormal      | +      | +            | -       |
| barcode03 | 3          | 51.6        | flox_deletion | -      | -            | -       |

#### Details.pdf

`Details.pdf`は上記CSVを可視化した以下のような図です.  

<img src="https://github.com/akikuno/DAJIN/blob/master/misc/images/Details.png" width="75%">  

barcode01は野生型コントロールです. barcode02とbarcode03はfloxノックインのゲノム編集を施したファウンダーマウスの結果です.  
barcode02のほぼ全てのアレルがintact target （flox以外の異常な変異の入っていないアレル）であることから、このマウスは目的のfloxアレルをホモでもつマウスの候補と考えられます.  

#### Consensus

`Conseusus` フォルダーには各アレルのコンセンサス配列が保存されています.  
ファイル形式はFASTAおよびHTMLです.  

HTMLでは色付けされた変異情報が表示されます.  

<a href="https://htmlpreview.github.io/?https://github.com/akikuno/DAJIN/blob/master/misc/images/tyr_c140cg.html" target= _blank rel= noopener> こちらは点変異のコンセンサス配列です. </a>

#### BAM

`BAM` フォルダーには解析したサンプルの全アレルおよび各アレルごとのBAMファイルが保存されています.  
この`BAM`ファイルは[IGV](http://software.broadinstitute.org/software/igv/)で可視化できます.  

