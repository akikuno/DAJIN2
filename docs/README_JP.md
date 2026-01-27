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
- メモリ推奨：16GB以上

>[!NOTE]
> 前身のDAJINは深層学習を使用していたため、効率的な計算にはGPUが必要でした。  
> 一方、**DAJIN2では深層学習を用いておらず、GPUは不要**です。  
> そのため、一般的なノートパソコンでも問題なく動作します。

### ソフトウェア

- Python 3.10-3.12
- Unix環境（Linux、macOS、WSL2など）

>[!IMPORTANT]
> **Windowsユーザーの方へ**  
> DAJIN2はLinux環境での実行を前提としています。Windowsをご利用の場合は、  
> WSL2（Windows Subsystem for Linux 2）を使用してください。

## [Bioconda](https://anaconda.org/bioconda/DAJIN2) （推奨）

```bash
# Setting up Bioconda
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
> ➡️ **最新バージョンの確認はこちら：**  
> https://github.com/akikuno/DAJIN2/releases
>
> 🔄 **最新バージョンへのアップデート方法：**
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

### 1. サンプルおよびコントロールファイル

DAJIN2では、ゲノム編集特異的な変異を検出するために、**ゲノム編集を受けていないコントロールサンプル**が必要です。  
ゲノム編集サンプルとコントロールのFASTQ/FASTA（gzip圧縮・非圧縮どちらも対応可能）、またはBAMファイルを含むディレクトリを指定します。

#### [Dorado](https://github.com/nanoporetech/dorado)によるベースコール

Doradoによるベースコール（[`dorado demux`](https://github.com/nanoporetech/dorado?tab=readme-ov-file#barcode-classification)）では、BAMファイルが出力されます：

```text
bam_pass
├── barcode01
│   └── EXP-PBC096_barcode01.bam
├── barcode02
│   └── EXP-PBC096_barcode02.bam
├── ...
└── unclassified
│   └── EXP-PBC096_unclassified.bam
```

> [!IMPORTANT]
> 各bamファイルを別々のディレクトリに格納してください。ディレクトリ名は任意です。


[`dorado correct`](https://github.com/nanoporetech/dorado#read-error-correction)によるシークエンスエラー補正後に出力されるFASTAファイルも同様に、別々のディレクトリに格納してください。

```text
dorado_correct
├── barcode01
│   └── EXP-PBC096_barcode01.fasta
└── barcode02
    └── EXP-PBC096_barcode02.fasta
```

> [!NOTE]
> Doradoによるベースコールでのマルチプレックスの使用方法については、[DORADO_HANDLING_JP.md](./DORADO_HANDLING_JP.md)をご覧ください。


#### [Guppy](https://community.nanoporetech.com/docs/prepare/library_prep_protocols/Guppy-protocol/v/gpb_2003_v1_revax_14dec2018/guppy-software-overview)によるベースコール

Guppyによるベースコール後は、以下のようなファイル構成で出力されます：

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


### 2. 想定アレル配列のFASTAファイル

FASTAファイルには、ゲノム編集によって想定されるアレルを記述します。

> [!IMPORTANT]
> コントロールアレルの指定：`>control`というヘッダー名とその配列は必須です。  

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

ここで、`>control`はコントロールアレルの配列を表しており、必須です。  
`>knock-in`と`>knock-out`はそれぞれノックインとノックアウトの想定アレル配列です。  

> [!IMPORTANT]
> **FASTA配列の両端は、アンプリコン配列の両端と一致させてください。**  
> アンプリコンより長い、または短い場合、その差分がIndelとして判定される可能性があります。  

## 単一サンプル解析

単一サンプルの解析コマンドは以下の通りです。  

```bash
DAJIN2 <-c|--control> <-s|--sample> <-a|--allele> <-n|--name> \
  [-g|--genome] [-b|--bed] [-t|--threads] [--no-filter] [-h|--help] [-v|--version]

Arguments:
  -c, --control             Path to the directory containing control FASTQ/FASTA/BAM files
  -s, --sample              Path to the directory containing sample FASTQ/FASTA/BAM files
  -a, --allele              Path to the FASTA file containing expected alleles
  -n, --name (optional)     Output directory name [default: Results]
  -g, --genome (optional)   Reference genome ID (e.g., hg38, mm39) [default: '']
  -b, --bed (optional)      Path to the BED6 file containing genome coordinates [default: '']
  -t, --threads (optional)  Number of threads to use [default: 1]
  --no-filter (optional)    Disable minor-allele filtering (keep alleles below 0.5%) [default: False]
  -h, --help                Show help message
  -v, --version             Show version
```

### 実行例

```bash
# Download the example dataset
curl -LJO https://github.com/akikuno/DAJIN2/raw/main/examples/example_single.tar.gz
tar -xf example_single.tar.gz

# Run DAJIN2 (single-sample analysis)
DAJIN2 \
    --control example_single/control \
    --sample example_single/sample \
    --allele example_single/stx2_deletion.fa \
    --name stx2_deletion \
    --bed example_single/stx2_deletion.fa \
    --threads 4
```


### BEDファイルを用いたゲノム座標の指定

参照ゲノムがUCSC提供のものではない場合、あるいはDAJIN2が依存する外部サーバー（UCSC Genome BrowserおよびGGGenome）が停止している場合は、`-b/--bed`オプションでBEDファイルを指定することで、オフラインで動作させられます。  

>[!IMPORTANT]
> UCSCやGGGenomeに依存する`--genome`よりも、安定して動作する`--bed`の利用を推奨します。

`-b/--bed`オプションを使用する際の注意点は以下の通りです：

**必ずBED6形式**を使用してください：

```
chr1    1000000    1001000    mm39    248956422    +
```

**各列の説明：**
- 1列目：染色体名（例：chr1、chr2）
- 2列目：開始位置（0インデックス）
- 3列目：終了位置（0インデックス）
- 4列目：名前（**ゲノムID**）
- 5列目：スコア（**染色体サイズを指定**）
- 6列目：鎖（`+` または `-`、**FASTAアレル配列の向きと一致させる必要あり**）

> [!NOTE]  
> 5列目のスコアには、1列目で指定した染色体のサイズを入力してください。  
> BED形式本来のスコアの上限は1000ですが、本ツールでは**染色体サイズを記載して問題ありません**。  

> [!IMPORTANT]  
> BEDファイルの6列目（鎖：`+`または`-`）は、**参照ゲノム配列に対するFASTAアレル配列の鎖の向きに合わせてください**。  
> - 参照配列が「ACCG」、FASTA配列が「ACCG」 → **フォワード鎖（5'→3'）**：`+`  
> - 参照配列が「ACCG」、FASTA配列が「CGGT」 → **リバース鎖（3'→5'）**：`-`

>[!NOTE]
> 詳細なBEDファイルの使用方法については、[BED_COORDINATE_USAGE.md](./BED_COORDINATE_USAGE.md)をご覧ください。

### `--no-filter`による希少変異の検出

DAJIN2は標準設定では、ノイズの軽減と精度向上のため、0.5%未満（100,000リードのダウンサンプリング中の5リード未満）のアレルを除外しています。  
しかし、希少変異や体細胞モザイクなど、非常に低頻度でマイナーアレルが存在する可能性がある場合は、`--no-filter`オプションを使用してこのフィルタリングを無効化できます。  

**`--no-filter`を使用する場面：**
- 希少体細胞変異の検出（<0.5%の頻度）
- 低レベルモザイクが疑われるサンプルの解析
- 頻度に関係なくすべての可能なアレルの検出が必要な研究

**使用方法：**
```bash
DAJIN2 \
    --control example_single/control \
    --sample example_single/sample \
    --allele example_single/stx2_deletion.fa \
    --name stx2_deletion \
    --bed example_single/stx2_deletion.fa \
    --threads 4 \
    --no-filter
```

> [!CAUTION]
> `--no-filter`の使用は、結果にノイズや偽陽性を増加させる可能性があります。希少アレルについては、追加の実験手法による検証を推奨します。


## 複数サンプルの一括処理

`batch`サブコマンドを利用することで、複数サンプルの同時処理が可能です。  
サンプル情報をまとめたCSVファイルまたはExcelファイルが必要です。  

> [!NOTE]
> サンプル情報のまとめ方は、[こちら](https://docs.google.com/presentation/d/e/2PACX-1vQMpqzwI9gtGnmMvqh9UFNxmpKDxcnUg74_TgLmd0FbBrrGQTa7CAQZvFlGDC2vxw/embed?start=false&loop=false&delayms=3000)をご参照ください。  

**必須列：** `sample`, `control`, `allele`, `name`  
**オプション列：** `genome`, `bed`（または`genome_coordinate`）、その他のカスタム列

**BEDファイルを使用したCSVの例：**
```csv
sample,control,allele,name,bed
/path/to/sample1,/path/to/control1,/path/to/allele1.fa,experiment1,/path/to/coords1.bed
/path/to/sample2,/path/to/control2,/path/to/allele2.fa,experiment2,/path/to/coords2.bed
```

> [!TIP]
> **同じ実験に属するサンプルには、`name`列に同じ値を使用することを推奨します。**  
> 同一の名前を使用することで、処理が並列化され、効率が向上します。  
> こちらが一例です 👉 [batch.csv](https://github.com/akikuno/DAJIN2/blob/main/examples/example_batch/batch.csv)

```bash
DAJIN2 batch <-f|--file> [-t|--threads] [--no-filter] [-h]

Arguments:
  -f, --file                Path to the CSV or Excel file
  -t, --threads (optional)  Number of threads to use [default: 1]
  --no-filter (optional)    Disable minor-allele filtering (keep alleles below 0.5%) [default: False]
  -h, --help                Show help message
```

### 実行例

```bash
# Download the example dataset
curl -LJO https://github.com/akikuno/DAJIN2/raw/main/examples/example_batch.tar.gz
tar -xf example_batch.tar.gz

# Run DAJIN2 (batch processing)
DAJIN2 batch --file example_batch/batch.csv --threads 4
```

## GUI（Graphical User Interface）モード

DAJIN2は、1コマンドで起動できるWebインターフェースを提供しています：

```bash
DAJIN2 gui
```

実行すると、既定のWebブラウザが開き、`http://localhost:{PORT}/`に以下のようなGUIが表示されます。

<img src="https://raw.githubusercontent.com/akikuno/DAJIN2/refs/heads/main/image/dajin2-gui.jpg" width="75%">

> [!NOTE]
> ブラウザが自動的に起動しない場合は、お使いのブラウザで`http://localhost:{PORT}/`に手動でアクセスしてください。


### GUIによる単一サンプル解析

1. **GUIの起動**  
   `DAJIN2 gui`を実行してWebインターフェースを開きます。

2. **プロジェクトの設定**  
   - **プロジェクト名**：任意の解析名を入力  
   - **ディレクトリアップロード**：サンプルまたはコントロールのFASTQ/FASTA/BAMファイルが含まれるディレクトリを選択  
   - **アレルFASTA**：想定アレル配列を含むFASTAファイルをアップロード  
   - **BEDファイル（任意）**：ゲノム座標を指定するBED6形式ファイルをアップロード

3. **パラメータの設定**  
   - **参照ゲノム（任意）**：UCSCゲノムID（例：`hg38`、`mm39`）を指定  
   - **スレッド数**：使用するCPUスレッド数を指定  
   - **フィルタ無効化**：0.5%未満の希少変異も検出する場合に有効化

4. **解析の実行**  
   「Start Analysis」をクリックすると、リアルタイムで進捗状況が表示されます。

5. **結果の確認**  
   解析完了後、出力フォルダのパスが表示され、結果ファイルにアクセスできます。


### GUIによるバッチ処理

1. **バッチファイルの準備**  
   `sample`、`control`、`allele`、`name` の列を含むCSVまたはExcelファイルを作成します。

2. **バッチファイルのアップロード**  
   「バッチ処理」タブから、作成した設定ファイルをアップロードします。

3. **全体設定の指定**  
   スレッド数やフィルタリングオプションを一括で設定します。

4. **進捗の監視**  
   各サンプルの解析状況が詳細なログ出力として表示されます。

5. **結果の確認**  
   結果は`DAJIN_Results/`フォルダ内に、サンプルごとのサブディレクトリとして保存されます。

# 📈 出力結果

DAJIN2の処理が完了すると、`DAJIN_Results/{NAME}`というディレクトリが作られます。  
この`DAJIN_Results/{NAME}`ディレクトリには、以下のファイルが含まれます：  

```text
DAJIN_Results/tyr-substitution
├── BAM
│   ├── control
│   ├── tyr_c230gt_01
│   ├── tyr_c230gt_10
│   └── tyr_c230gt_50
├── DAJIN2_log_20260127_140954_076887.txt
├── FASTA
│   ├── tyr_c230gt_01
│   ├── tyr_c230gt_10
│   └── tyr_c230gt_50
├── HTML
│   ├── tyr_c230gt_01
│   ├── tyr_c230gt_10
│   └── tyr_c230gt_50
├── MUTATION_INFO
│   ├── tyr_c230gt_01.csv
│   ├── tyr_c230gt_10.csv
│   └── tyr_c230gt_50.csv
├── VCF
│   ├── tyr_c230gt_01
│   ├── tyr_c230gt_10
│   └── tyr_c230gt_50
├── launch_report_mac.command
├── launch_report_windows.bat
└── read_summary.xlsx
```

## 1. launch_report_windows.bat / launch_report_mac.command

Windowsの場合は`launch_report_windows.bat`、macOSの場合は`launch_report_mac.command`をダブルクリックすると、  
ブラウザが起動し、レポートが表示されます。

<img src="https://raw.githubusercontent.com/akikuno/logos/refs/heads/main/DAJIN2/DAJIN2-report.jpg" width="100%" />

>[!IMPORTANT]
> **積み上げ棒グラフにて興味のあるアレルをクリックすると、配列にどのような変異が生じているのかの詳細を閲覧できます（上図右）。**

図中の**Allele type**はアレルの種類を、**Percent of reads**は該当するリードのアレル割合を示しています。  

**Allele type**の種類は以下の通りです：  
- **{アレル名}**：ユーザーがFASTAにより定義したアレルに完全一致するアレル
- **{アレル名} with indels**：ユーザーがFASTAにより定義したアレルに類似しているが、数塩基の置換、欠失、挿入、逆位を含むアレル
- **unassigned insertion/deletion/inversion**：ユーザーが定義していない10塩基以上の欠失、挿入、逆位を含むアレル



> [!WARNING]  
> PCRアンプリコンを用いたターゲットシーケンシングでは、増幅バイアスのため**Percent of reads**がサンプルに存在するアレルの割合と一致しないことがあります。  
> とくに大型欠失が存在する場合、欠失アレルが顕著に増幅されることから、実際のアレル割合を反映しない可能性が高まります。  

## 2. read_summary.xlsx

read_summary.xlsxには、各アレルのリード数と存在割合が記載されています。  
レポートの積み上げ棒グラフは、`read_summary.xlsx`を可視化したものです。  
論文用の図を作成する際などの資料としてご利用ください。  


## 3. BAMとVCF

BAMおよびVCFディレクトリには、アレルごとに分類されたBAMファイルおよびVCFファイルが格納されています。  

> [!NOTE]  
> `--bed`や`--genome`による指定がない場合、入力のFASTAファイルのcontrolアレルにアライメントされます。

## 4. FASTAとHTML

FASTAディレクトリには、各アレルのFASTAファイルが保存されます。  
HTMLディレクトリには、変異箇所が色付けされた各アレルのHTMLファイルが保存されます。  
Tyr点変異（緑色）の例を以下に示します：

<img src="https://raw.githubusercontent.com/akikuno/logos/refs/heads/main/DAJIN2/tyr-substitution.png" width="75%" />

また、DAJIN2はサンプルに含まれる代表的なSVアレル（Insertion, Deletion, Inversion）を抽出し、SV箇所を下線で色付けします。  
以下は、Inversion（下線紫色）の両端に欠失（水色）および挿入（赤色）が認められる例です。  

<img src="https://raw.githubusercontent.com/akikuno/logos/refs/heads/main/DAJIN2/cables2-inversion.png" width="75%" />


## 5. MUTATION_INFO

MUTATION_INFOディレクトリには、各アレルの変異箇所を示すテーブルが保存されます。  
*Tyr*点変異の例を以下に示します：
- 点変異の染色体上の位置と、変異の種類が記載されています。

<img src="https://user-images.githubusercontent.com/15861316/274519342-a613490d-5dbb-4a27-a2cf-bca0686b30f0.png" width="75%">


# 📣 お問い合わせ

ご質問やバグ報告、ご意見・ご要望など、皆さまのフィードバックをお待ちしています。  

ご報告は以下のGoogleフォームからお願いいたします。  
👉 [Googleフォーム](https://forms.gle/r4YRs1th7NGHfDcS9)

また、GitHubアカウントをお持ちの方は、Issueからもご報告いただけます（日本語可）。  
👉 [GitHub Issues](https://github.com/akikuno/DAJIN2/issues/new/choose)  

> [!NOTE]
> よくあるご質問については、[こちら](https://github.com/akikuno/DAJIN2/blob/main/docs/FAQ_JP.md)をご覧ください。


# 🤝 行動規範

本プロジェクトは[Contributor Code of Conduct（コントリビューター行動規範）](https://github.com/akikuno/DAJIN2/blob/main/docs/CODE_OF_CONDUCT.md)に基づいて公開されています。  

# 📄 参考文献


[Kuno A, et al. (2022) DAJIN enables multiplex genotyping to simultaneously validate intended and unintended target genome editing outcomes. *PLoS Biology* 20(1): e3001507.](https://doi.org/10.1371/journal.pbio.3001507)
