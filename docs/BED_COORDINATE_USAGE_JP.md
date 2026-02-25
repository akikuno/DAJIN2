# BEDファイルを用いたゲノム座標の指定

このドキュメントでは、`-b/--bed`オプションで必要となるBEDファイルについて解説します。  

---

## 概要

DAJIN2は現在、`-b/--bed`オプション（後方互換性のため`--genome-coordinate`も可）を用いて、BEDファイルによるゲノム座標指定をサポートしています。これにより、外部サーバー（UCSC Genome BrowserおよびGGGENOME）が利用できない場合でも、オフラインで動作可能になります。

> [!IMPORTANT]
> `--genome`オプションも有効ですが、外部サーバーとの通信エラーの可能性があるため、`-b/--bed`のご利用を推奨します。


### 優先ルール

1. `-b/--bed`が指定されている場合、`--genome`より**優先**されます
2. 両方が指定されている場合、BEDファイルの座標が使用されます
3. `--genome`のみが指定されている場合、従来どおり外部サーバー（UCSC Genome BrowserおよびGGGENOME）との通信を介してゲノム座標を入手します
4. どちらも指定されていない場合、ゲノム座標なしで解析が実行されます

### BEDファイル形式

DAJIN2は、以下の要件を満たすBED6形式ファイルを受け付けます。

````
chr1	1000000	1001000	mm39	195154279	+
````

#### カラム定義

- **1列目**: 染色体名（例:`chr1`,`chr2`,`chrX`または`1`,`2`,`X`）
- **2列目**: 開始位置（0インデックス）
- **3列目**: 終了位置（0インデックス）
- **4列目**: 名前（**ゲノムID**、例:`mm39`,`hg38`）
- **5列目**: スコア（**IGVの適切な可視化のための染色体サイズ**、例:chr1の`248956422`）
- **6列目**: ストランド（`+`または`-`、**FASTAアレル配向と一致している必要があります**）

> [!NOTE]  
> scoreフィールド（5列目）には、1列目で指定した染色体のサイズを入力してください。    
> オリジナルのBED形式ではscoreは1000までに制限されていますが、DAJIN2では**染色体サイズを問題なく受け付けます**。  

> [!IMPORTANT]  
> **ストランド方向は一致している必要があります**。BEDファイルのストランドフィールド（6列目:`+`または`-`）は、**FASTAアレル配列のストランド方向と一致している必要があります**。  
> - FASTAアレル配列が**順鎖（forward strand）**（5'→3'）である場合、BEDファイルでは`+`を使用してください  
> - FASTAアレル配列が**逆鎖（reverse strand）**（3'→5'）である場合、BEDファイルでは`-`を使用してください


**染色体サイズの調べ方:**

- Human（hg38）: https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&chromInfoPage=
- Mouse（mm39）: https://genome.ucsc.edu/cgi-bin/hgTracks?db=mm39&chromInfoPage=
- または:`samtools faidx genome.fa`を使用してFASTAインデックスから染色体サイズを取得

#### ストランド処理

- **ストランドフィールド（6列目）** は、適切な配列方向のために必須です
- **`+`ストランド**: 配列は順鎖（5'→3'）上にあります
- **`-`ストランド**: 配列は逆鎖上にあり、逆相補鎖に変換されます（reverse complement）
- **必須形式**: DAJIN2はストランド情報（`+`または`-`）を含むBED6形式を必要とします
- **IGV可視化**: ストランド情報はBAMファイル内に保持され、ゲノムブラウザ表示に利用されます
- **配列処理**: DAJIN2はマイナスストランド領域に対して自動的にreverse complementを適用します

### ファイル拡張子

サポートされるファイル拡張子:

- `.bed`
- `.bed.gz`

### 複数区間を含むBEDファイル

- BEDファイルに複数の区間が含まれている場合、DAJIN2は解析に**1行目の区間**を使用します
- 複数区間が検出された場合は警告ログが出力されます
- 将来のバージョンでは複数領域解析をサポートする可能性があります

---

## 使い方

### 単サンプル解析

````bash
DAJIN2 --sample sample/ \
    --control control/ \
    --allele alleles.fa \
    --name experiment1 \
    --bed coordinates.bed
````

### バッチモード対応

`-b/--bed`オプションはバッチモードでもサポートされています。

### バッチファイル形式

````csv
sample,control,allele,name,genome,bed
sample1/,control/,alleles.fa,exp1,hg38,exp1.bed
sample2/,control/,alleles.fa,exp2,,exp2.bed
````

- 各行で`bed`列に独自のBEDファイルを指定できます
- 同じ優先ルールが適用されます:`bed`は`genome`より優先されます


---

## エラーハンドリング

よくあるエラーと解決策:

**エラー**: `BED file must have .bed or .bed.gz extension`

- **解決策**: 
- 
- ファイル名を正しい拡張子に変更してください

**エラー**: `DAJIN2 requires BED6 format with strand information`

- **解決策**: BEDファイルが6列あり、6列目にストランド（`+/-`）を含んでいることを確認してください

**エラー**: `Invalid chromosome size format at line X: 'Y' (must be an integer)`

- **解決策**: 5列目を染色体サイズを表す正の整数に設定してください  

**エラー**: `Invalid chromosome size at line X: Y (must be a positive integer)`

- **解決策**: 5列目に0より大きい正の整数が入っていることを確認してください

**エラー**: `Invalid or missing strand at line X: 'Y' (must be '+' or '-')`

- **解決策**: 6列目をストランド情報として`+`または`-`のいずれかに設定してください

**エラー**: `Invalid end position at line X: Y (must be > start)`

- **解決策**: 終了位置が開始位置より大きいことを確認してください

**エラー**: `No valid intervals found in BED file`

- **解決策**: BEDファイルに少なくとも1つの有効なゲノム座標が含まれていることを確認してください


---

## 後方互換性

- 既存の`--genome`機能はすべて維持されます
- `genome_coordinate`列を持たない既存のバッチファイルも引き続き動作します
- 出力形式やファイル構造への変更はありません

