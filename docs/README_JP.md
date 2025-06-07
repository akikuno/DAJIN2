[![License](https://img.shields.io/badge/License-MIT-9cf.svg)](https://choosealicense.com/licenses/mit/)
[![Test](https://img.shields.io/github/actions/workflow/status/akikuno/dajin2/pytest.yml?branch=main&label=Test&color=brightgreen)](https://github.com/akikuno/dajin2/actions)
[![Python](https://img.shields.io/pypi/pyversions/DAJIN2.svg?label=Python&color=blue)](https://pypi.org/project/DAJIN2/)
[![PyPI](https://img.shields.io/pypi/v/DAJIN2.svg?label=PyPI&color=orange)](https://pypi.org/project/DAJIN2/)
[![Bioconda](https://img.shields.io/conda/v/bioconda/dajin2?label=Bioconda&color=orange)](https://anaconda.org/bioconda/dajin2)
[![DOI](https://zenodo.org/badge/387721337.svg)](https://zenodo.org/badge/latestdoi/387721337)
[![Paper](https://img.shields.io/badge/Plos%20Biol-10.1371/journal.pbio.3001507-lightgreen)](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3001507)
[![お問い合わせ](https://img.shields.io/badge/お問い合わせ-923DE2)](https://forms.gle/r4YRs1th7NGHfDcS9)


<p align="center">
<img src="https://user-images.githubusercontent.com/15861316/261833016-7f356960-88cf-4574-87e2-36162b174340.png" width="90%">
</p>

DAJIN2は、ナノポアターゲットシーケンシングを用いた、ゲノム編集サンプルの遺伝型解析ツールです。

**DAJIN**の名称は「一網**打尽**」に由来しており、意図された変異だけでなく、意図しない変異も含めて、**ゲノム編集の結果を一括かつ網羅的に検出する**という本ツールの設計思想が込められています。


# 🌟 特徴

+ **網羅的な変異検出**  
  ナノポアによるターゲット領域の解析を通じて、点変異から構造多型まで、あらゆるゲノム編集イベントを検出可能です。  
  特に、**想定外の変異**や、**欠失領域における挿入**など、複合的な変異の検出に優れています。

+ **高感度なアレル分類**  
  モザイクアレルの分類に対応しており、約1%のマイナーアレルも検出可能です。

+ **直感的な可視化**  
  ゲノム編集の結果を視覚的に表示し、変異の迅速かつ容易な識別を支援します。

+ **多サンプル対応**  
  複数サンプルの一括処理に対応しており、大規模な実験や比較研究を効率的に実施できます。

+ **簡単なインストールと操作**  
  特別な計算機環境を必要とせず、一般的なノートパソコンでも動作します。  
  BiocondaやPyPIから簡単に導入でき、コマンドラインから直感的に利用可能です。  


# 🛠 インストール

## 環境

### ハードウェア

- **一般的なノートパソコンで動作可能**
- メモリ推奨：8GB以上

>[!NOTE]
> 前身のDAJINは深層学習を使用していたため、効率的な計算にはGPUが必要でした。  
> 一方、**DAJIN2では深層学習を用いておらず、GPUは不要**です。  
> そのため、一般的なノートパソコンでも問題なく動作します。

### ソフトウェア

- Python >= 3.9
- Unix環境 (Linux, macOS, WSL2, etc.)

>[!IMPORTANT]
> **Windowsユーザーの方へ**  
> DAJIN2はLinux環境での実行を前提としています。Windowsをご利用の場合は、  
> **WSL2（Windows Subsystem for Linux 2）** を使用してください。

## [Bioconda](https://anaconda.org/bioconda/DAJIN2) （推奨）

```bash
# Setup of Bioconda
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority flexible

# Install DAJIN2
conda create -n env-dajin2 python=3.12 DAJIN2 -y
conda activate env-dajin2
```

## [PyPI](https://pypi.org/project/DAJIN2/)

```bash
pip install DAJIN2
```

> [!IMPORTANT]
> **DAJIN2は活発に開発・改良が進められています。**  
> 最新機能をご利用いただくために、ご使用中のバージョンが最新版であることをご確認ください。  
>
> 🔍 **現在のバージョンの確認方法：**
> ```bash
> DAJIN2 --version
> ```
>
> 🆕 **最新版の確認はこちら：**  
> https://github.com/akikuno/DAJIN2/releases
>
> 🔄 **最新版へのアップデート方法：**
> ```bash
> conda update DAJIN2 -y
> ```
> または
> ```bash
> pip install -U DAJIN2
> ```


> [!CAUTION]
> インストールに問題が発生した場合は、[トラブルシューティングガイド](https://github.com/akikuno/DAJIN2/blob/main/docs/TROUBLESHOOTING.md)をご覧ください。



# 💻 使用方法

## 必要なファイル

### サンプルおよびコントロールファイル

DAJIN2では、ゲノム編集特異的な変異を検出するために、**ゲノム編集を受けていないコントロールサンプル**が必要です。  
ゲノム編集サンプルとコントロールのFASTQ/FASTA（gzip圧縮・非圧縮どちらも対応可能）、またはBAMファイルを含むディレクトリを指定します。

#### [Guppy](https://community.nanoporetech.com/docs/prepare/library_prep_protocols/Guppy-protocol/v/gpb_2003_v1_revax_14dec2018/guppy-software-overview)によるベースコール

Guppyによるベースコール後、以下のようなファイル構成が出力されます：

```text
fastq_pass
├── barcode01
│   ├── fastq_runid_b347657c88dced2d15bf90ee6a1112a3ae91c1af_0_0.fastq.gz
│   ├── fastq_runid_b347657c88dced2d15bf90ee6a1112a3ae91c1af_10_0.fastq.gz
│   └── fastq_runid_b347657c88dced2d15bf90ee6a1112a3ae91c1af_11_0.fastq.gz
└── barcode02
    ├── fastq_runid_b347657c88dced2d15bf90ee6a1112a3ae91c1af_0_0.fastq.gz
    ├── fastq_runid_b347657c88dced2d15bf90ee6a1112a3ae91c1af_10_0.fastq.gz
    └── fastq_runid_b347657c88dced2d15bf90ee6a1112a3ae91c1af_11_0.fastq.gz
```

上記のbarcode01をコントロール、barcode02をサンプルとすると、それぞれのディレクトリは下記の通りに指定します：

+ コントロール: `fastq_pass/barcode01`
+ サンプル: `fastq_pass/barcode02`

#### [Dorado](https://github.com/nanoporetech/dorado)によるベースコール

Doradoによるベースコール（[`dorado demux`](https://github.com/nanoporetech/dorado?tab=readme-ov-file#barcode-classification)）においては、BAMファイルが出力されます：

```text
dorado_demultiplex
├── EXP-PBC096_barcode01.bam
└── EXP-PBC096_barcode02.bam
```

> [!IMPORTANT]
> 各bamファイルを別々のディレクトリに格納してください。ディレクトリ名は任意です。

```text
dorado_demultiplex
├── barcode01
│   └── EXP-PBC096_barcode01.bam
└── barcode02
    └── EXP-PBC096_barcode02.bam
```

[`dorado correct`](https://github.com/nanoporetech/dorado)によるシークエンスエラー補正後に出力されるFASTAファイルも同様に、別々のディレクトリに格納してください。

```text
dorado_correct
├── barcode01
│   └── EXP-PBC096_barcode01.fasta
└── barcode02
    └── EXP-PBC096_barcode02.fasta
```

上記のbarcode01をコントロール、barcode02をサンプルとすると、それぞれのディレクトリは下記の通りに指定します：

+ コントロール: `dorado_demultiplex/barcode01` / `dorado_correct/barcode01`
+ サンプル: `dorado_demultiplex/barcode02` / `dorado_correct/barcode02`


### 想定アレル配列のFASTAファイル

FASTAファイルには、ゲノム編集によって想定されるアレルを記述します。

> [!IMPORTANT]
>コントロールアレルの指定： `>control`というヘッダー名と、その配列は必須です。  

事前に想定されるアレルがある場合（例：ノックインやノックアウト）、それらのシーケンスもFASTAファイルに記載してください。これらの想定アレルの名称は任意に設定できます。

以下は、FASTAファイルの一例です：

```text
>control
ACGTACGTACGTACGT
>knock-in
ACGTACGTCCCCACGTACGT
>knock-out
ACGTACGT
```

ここで、`>control` はコントロールアレルの配列を表しており、必須です。  
`>knock-in` と `>knock-out` はそれぞれノックインとノックアウトの想定アレル配列です。  

> [!IMPORTANT]
> **FASTA配列の両端は、アンプリコン配列の両端と一致させてください。**   
> アンプリコンよりも長い、または短い場合、その差分はIndelとして判定される可能性があります  

## 単一サンプル解析

単一サンプルの解析コマンドは以下の通りです。

```bash
DAJIN2 <-c|--control> <-s|--sample> <-a|--allele> <-n|--name> \
  [-g|--genome] [-t|--threads] [-h|--help] [-v|--version]

引数:
  -c, --control             コントロールのFASTQ/FASTA/BAMファイルが格納されたディレクトリのパス
  -s, --sample              サンプルのFASTQ/FASTA/BAMファイルが格納されたディレクトリのパス
  -a, --allele              ゲノム編集によって想定されるアレルを記載したFASTAファイルのパス
  -n, --name (オプション)     出力ディレクトリの名前 [デフォルト: Results]
  -g, --genome (オプション)   参照ゲノムのID (e.g hg38, mm39) [デフォルト: '']
  -t, --threads (オプション)  使用するスレッド数 [デフォルト: 1]
  -h, --help                ヘルプメッセージの出力
  -v, --version             バージョンの出力
```

### 実行例

```bash
# Example datasetのダウンロード
curl -LJO https://github.com/akikuno/DAJIN2/raw/main/examples/example_single.tar.gz
tar -xf example_single.tar.gz

# DAJIN2の実行（単一サンプル解析）
DAJIN2 \
    --control example_single/control \
    --sample example_single/sample \
    --allele example_single/stx2_deletion.fa \
    --name stx2_deletion \
    --genome mm39 \
    --threads 4
```

## 複数サンプルの一括処理

`batch`サブコマンドを利用することで、複数サンプルの同時処理が可能です。  
サンプル情報をまとめたCSVファイルまたはExcelファイルが必要です。  

> [!NOTE]
> サンプル情報のまとめ方は、[こちら](https://docs.google.com/presentation/d/e/2PACX-1vQMpqzwI9gtGnmMvqh9UFNxmpKDxcnUg74_TgLmd0FbBrrGQTa7CAQZvFlGDC2vxw/embed?start=false&loop=false&delayms=3000)をご参照ください。  

> [!TIP]
> **同じ実験に属するサンプルには、`name`列に同じ値を使用することを推奨します。**  
> 同一の名前を使用することで、処理が並列化され、効率が向上します。  
> こちらが一例です 👉 [batch.csv](https://github.com/akikuno/DAJIN2/blob/main/examples/example_batch/batch.csv)

```bash
DAJIN2 batch <-f|--file> [-t|--threads] [-h]

引数:
  -f, --file                CSVまたはExcelファイルのパス
  -t, --threads (オプション)  使用するスレッド数 [デフォルト: 1]
  -h, --help                ヘルプメッセージの出力
```

### 実行例

```bash
# Example datasetのダウンロード
curl -LJO https://github.com/akikuno/DAJIN2/raw/main/examples/example_batch.tar.gz
tar -xf example_batch.tar.gz

# DAJIN2の実行（バッチ処理）
DAJIN2 --file batch.csv --threads 4
```


# 📈 出力結果

DAJIN2の処理が完了すると、`DAJIN_Results/{NAME}`というディレクトリが作られます。  
この`DAJIN_Results/{NAME}`ディレクトリには、以下のファイルが含まれます：  

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
├── read_plot.html
├── read_plot.pdf
└── read_summary.xlsx
```

## 1. BAM

BAMディレクトリには、アレルごとに分類されたBAMファイルが格納されています。  

> [!NOTE]  
> `genome`オプションで参照ゲノムを指定すると、その参照ゲノムにアライメントされます。  
> 指定がない場合、入力のFASTAファイルのcontrolアレルにアライメントされます。

## 2. FASTA と HTML

FASTAディレクトリには、各アレルのFASTAファイルが保存されます。  
HTMLディレクトリには、変異箇所が色付けされた各アレルのHTMLファイルが保存されます。  
Tyr点変異（緑色）の例を以下に示します：


<img src="https://raw.githubusercontent.com/akikuno/logos/refs/heads/main/DAJIN2/tyr-substitution.png" width="75%" />


また、DAJIN2はサンプルに含まれる代表的なSVアレル(Insertion, Deletion, Inversion)を抽出し、SV箇所を下線で色付けします。  
以下は、Inversion（下線紫色）の両端に欠失（水色）および挿入（赤色）が認められる例です：  

<img src="https://raw.githubusercontent.com/akikuno/logos/refs/heads/main/DAJIN2/cables2-inversion.png" width="75%" />


## 3. MUTATION_INFO

MUTATION_INFOディレクトリには、各アレルの変異箇所を示すテーブルが保存されます。  
*Tyr*点変異の例を以下に示します：
- 点変異の染色体上の位置と、変異の種類が記載されています。

<img src="https://user-images.githubusercontent.com/15861316/274519342-a613490d-5dbb-4a27-a2cf-bca0686b30f0.png" width="75%">

## 4. read_summary.xlsx / read_plot.html / read_plot.pdf

read_summary.xlsxには、各アレルのリード数と存在割合が記述されています。  
read_plot.html および read_plot.pdf は、resd_summary.xlsxを可視化したもので、各アレルの割合を図示しています。  
図中の**Allele type**はアレルの種類を、**Percent of reads**は該当するリードのアレル割合を示しています。  

**Allele type**の種類は以下の通りです：

- **Intact**：入力のFASTAアレルと完全に一致するアレル
- **Indels**：50塩基以内の置換、欠失、挿入、逆位を含むアレル
- **SV**：50塩基以上の置換、欠失、挿入、逆位を含むアレル


<img src="https://user-images.githubusercontent.com/15861316/274521067-4d217251-4c62-4dc9-9c05-7f5377dd3025.png" width="75%">

> [!WARNING]  
> PCRアンプリコンを用いたターゲットシーケンシングでは、増幅バイアスのため **% of reads**が実際のアレルの割合と一致しないことがあります。  
> とくに大型欠失が存在する場合、欠失アレルが顕著に増幅されることから、実際のアレル割合を反映しない可能性が高まります。


# 📣 お問い合わせ

ご質問やバグ報告、ご意見・ご要望など、皆さまのフィードバックをお待ちしています。  

ご報告は、以下のGoogleフォームからお願いいたします：  
👉 [Googleフォーム](https://forms.gle/r4YRs1th7NGHfDcS9)

また、GitHubアカウントをお持ちの方は、Issueからもご報告頂けます（日本語可）：  
👉 [GitHub Issues](https://github.com/akikuno/DAJIN2/issues/new/choose)  

> [!NOTE]
> よくあるご質問については、[こちら](https://github.com/akikuno/DAJIN2/blob/main/docs/FAQ_JP.md)をご覧ください。


# 🤝 行動規範

本プロジェクトは [Contributor Code of Conduct（コントリビューター行動規範）](https://github.com/akikuno/DAJIN2/blob/main/docs/CODE_OF_CONDUCT.md)に基づいて公開されています。  

# 📄 参考文献


[Kuno A, et al. (2022) DAJIN enables multiplex genotyping to simultaneously validate intended and unintended target genome editing outcomes. *PLoS Biology* 20(1): e3001507.](https://doi.org/10.1371/journal.pbio.3001507)
