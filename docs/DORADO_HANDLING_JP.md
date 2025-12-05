# DoradoによるBase-callおよびDemultiplexの実行手順

> [!CAUTION]
> 以下の手順は、Nanopore社製のDoradoを用いたBase-callおよびDemultiplexの実行手順です。  
> 本ドキュメントは2025年12月時点の情報に基づいて作成されています。
> 最新のドキュメントは以下のURLを参照してください。  
> https://software-docs.nanoporetech.com/dorado/latest/  

## (1) Doradoのインストール

Doradoのバイナリをダウンロードし、実行権限を付与します。

```bash
curl -LO https://cdn.oxfordnanoportal.com/software/analysis/dorado-1.3.0-linux-x64.tar.gz
tar -xzf dorado-1.3.0-linux-x64.tar.gz
chmod +x ./dorado-1.3.0-linux-x64/bin/dorado
```

## (2) Base-callの実行

Doradoを用いてBase-callを実行します。以下のコマンドでは、`dna_r10.4.1_e8.2_400bps_sup`モデルを使用しています。適宜モデルを変更してください。

```bash
./dorado-1.3.0-linux-x64/bin/dorado basecaller \
    sup \
    pod5/ \
    --device cuda:all \
    --kit-name EXP-PBC096 \
    > calls.bam
```

> [!IMPORTANT]
> `--kit-name`オプションは、後述するDemultiplexの際に必要となります。  
> 適切なキット名を指定してください。  
> DAJIN論文では、`EXP-PBC096`を使用しています。

## (3) Demultiplexの実行

Doradoを用いてDemultiplexを実行します。以下のコマンドでは、Base-callで生成された`calls.bam`を入力としています。

```bash
mkdir -p dorado-demux

./dorado-1.3.0-linux-x64/bin/dorado demux \
    --output-dir dorado-demux \
    --no-classify \
    --threads 8 \
    calls.bam

```

## (4) Demultiplex結果の確認

Demultiplexの結果は、指定した出力ディレクトリ内に保存されます。  
各バーコードごとにフォルダが作成され、その中にBase-callされたBAMファイルが格納されています。

以下は、出力ディレクトリの例です。  

```bash
dorado-demux
└── 12635
    └── 20251110_0508_0_PBI23287_73af7da2
        └── bam_pass
            ├── barcode01
            │   └── PBI23287_pass_barcode01_73af7da2_00000000_0.bam
            ├── barcode02
            │   └── PBI23287_pass_barcode02_73af7da2_00000000_0.bam
            ├── barcode03
            │   └── PBI23287_pass_barcode03_73af7da2_00000000_0.bam
            └── unclassified
                └── PBI23287_pass_unclassified_73af7da2_00000000_0.bam
```

DAJIN2では、barcodeのディレクトリを入力ファイルとして利用します。  
例えば、batch.csvは以下のようになります。

| name | control                                                                        | sample                                                                         | allele  | bed      |
| ---- | ------------------------------------------------------------------------------ | ------------------------------------------------------------------------------ | ------- | -------- |
| test | dorado-demux-100000/12635/20251110_0508_0_PBI23287_73af7da2/bam_pass/barcode01 | dorado-demux-100000/12635/20251110_0508_0_PBI23287_73af7da2/bam_pass/barcode02 | test.fa | test.bed |
| test | dorado-demux-100000/12635/20251110_0508_0_PBI23287_73af7da2/bam_pass/barcode01 | dorado-demux-100000/12635/20251110_0508_0_PBI23287_73af7da2/bam_pass/barcode03 | test.fa | test.bed |

パスが長いので、`bam_pass`ディレクトリに移動してからDAJIN2を実行することをお勧めします。  
その際に必要なbatch.csvの例を以下に示します。  

| name | control   | sample    | allele  | bed      |
| ---- | --------- | --------- | ------- | -------- |
| test | barcode01 | barcode02 | test.fa | test.bed |
| test | barcode01 | barcode03 | test.fa | test.bed |

> [!NOTE]
> ディレクトリを移動した場合、`allele`および`bed`のファイルパスも適宜変更してください。
