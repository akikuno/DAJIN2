# ToDO

+ 🔴: Better to focus on.
+ 🟡: One day experiments.
+ 🟢: Need a time to complete.

## 開発一般

+ [ ] 🟡Controlを先に処理する。
  + `main`関数を`check_inputs`、`execute_control`、`execute_sample`の３つに分ける
  + `NAME`が同じで`CONTROL`が違うという場合も考慮する

+ [ ] 🟡debugを入れる
+ [ ] 🟡低スペックPCでメモリがオーバーしないか確認する
+ [ ] 🟢README.md, READMD_jp.mdの執筆
+ [ ] 🟢品質の低い（リード長が短く、バラバラとしている）Controlでもジェノタイピングができるか検討する
+ [ ] 🟢並列処理を組み入れる

## DAJIN2/core/main.py

###  前処理

+ [ ] 🟡`wslPath`をcore.mainのなかに入れる

### マッピング

### クラスタリング

+ [ ] 🔴Opticsでは非常に細かく分かれすぎる問題がある (AyabeTask1. barcode31など)
  + CRISPRESSO2のように細かく分かれていても見やすいような図示にする？
    + Plotlyに良さそうな関数がありました：https://plotly.com/python/alignment-chart/
  + 一度細かくクラスタリングしたあと、コントロールとの引き算を行う？
  + 共通領域を取り出してマージする？
+ [ ] (test) `clustering`

### コンセンサス

+ [ ] (test) call_percentage
+ [ ] (test) join_listdicts

### レポート

+ [ ] 🟢batchでアレル頻度の作図関数を作る
  + アレル割合が多い順に積み重ねたい
  + https://stackoverflow.com/questions/63920999/plotly-sorting-the-y-axis-bars-of-a-stacked-bar-chart-by-value

+ [ ] HTMLファイル
  + [ ] 🟢CRISPResso2のように割合もわかるとヘテロのときに嬉しい
  + [ ] 🟢NとDが重なった時にNが青色になるので、Nを優先させる

+ [ ] 🔴VCFファイル

+ [ ] 🟡BAMファイル
  + [ ] Hrでプラス鎖の確認
  + [ ] (test) `report_bam`
  + [x] controlも含める
  + [x] Tyrでマイナス鎖からプラス鎖に変換

+ [ ] 🟢名前をどうするか…
  + 2022-10-19時点では`barcode31_allele1_flox_intact_52.979%`のような感じになっている。
    + 視認性が良くない。
    + とくに`flox_intact`や`control_mutated`などのところ
    + `controlアレルに類似していて、かつcontrolアレルにはない変異が入っている`という意味をいかにして伝えるか
  + サンプル名をフォルダ名にする？


## DAJIN view

+ [ ] 🔴各サンプルごとにviewを分けたほうが見やすい
+ [ ] 🟢バックグラウンド実行
+ [ ] (test) `view.execute()`

## DAJIN gui

+ [ ] 🔴ローディングページ
+ [ ] 🔴終了時のインストラクションページ
+ [ ] 🟢バックグラウンド実行
+ [ ] 🟢デザインおよびレイアウト


## その他

+ [ ] 🟢サンプルのみ（コントロールなし）のときにも（可能な限り）アレルを分離できるようにしたい
+ [ ] 🟢`CSSPLIT`は名が体を表し切れていないので改名したい
+ [ ] 🟢ターゲットアンプリコンではないとき（染色体にマッピングするとき）、MIDSV変換が巨大な文字列となってしまうときの対応

---

# MEMO

## tyr

+ アルビノ点変異は829番目の塩基

## AyabeTask1

+ barcode33の1塩基欠失は`1748`番目の塩基

+ キメラDNAのleft-loxpは`control, SV false`にある。例: 880d5333-7a52-4dad-95fb-b778f06c1e5f
  + left-loxpがInsertion (at 1735) として扱われている。 left-loxpは46塩基長（<50）のためSVではない。
+ キメラDNAのright-loxpは`flox, SV false`にある。例：6f536ce8-32d6-41d5-8b09-66bef425b9d8
  + left-loxpがDeletionとして扱われる。 left-loxpは46塩基長（<50）のためSVではない。


---
# DONE

## 開発一般

+ [x] (2022-10-29) 点変異のヘテロ（ポジコン）と10%, 5%, 1%のサンプルを用意する（in vivoゲノム編集を想定）

+ [x] 引数チェック関数 → メイン関数 → レポート関数の3部構成とする
  + メイン関数はbatchでも使えるようにする
  + 引数チェックとレポートはbatchとsingleで異なる（可能性がある）
  + -> (2022-10-07) batch.py, single.pyのなかに引数チェック関数とレポート関数を入れる方式にする

+ [x] (2022-10-07) `.tempdir-{name}`ではなく`.tempdir/{name}`とする (2022-10-08)
+ [x] batchモードの搭載
+ [x] GitHub Actionsによるテスト自動化

+ [x] `is_control`の実装
  + batchモードの時にcontrolが計算済みであれば処理を省略する
  + controlのファイルサイズが同じで、midsvのjsonlがあればコントロールは処理済みとする
  + メイン関数内のcontrolとsampleの処理を分ける

## 前処理

+ [x] `DICT_ALLELE`は80文字改行のFASTAに対応しているのか確認する -> `mappy.fastx_read`で必ず1つにまとまることを確認した
+ [x] (test) `format_inputs`


+ [x] 入力ファイルのフォーマット (format_input.py)
  + [x] sample, controlのテキスト化
  + [x] alleleでmulti-fastaをsingle fastaにする
  + [x] alleleでwtかcontrolがない場合にはエラー終了
  + [x] alleleで重複配列がある場合にはエラー終了

+ [x] `check_inputs`と`format_inputs`を分ける
+ [x] (test) `check_inputs`
+ [x] 大型欠失とInversionは断端のマイクロホモロジーがあるので、配列とQualityを元にして結合する -> `midsv v0.7.0`で実装

## マッピング
+ [x] mappyによるマッピング

## クラスタリング

+ [x] コントロールとの引き算を検討する
  +  いまはコントロールと10%以上の差があるとDIFFLOCIとして出力されるが、それではin vivoゲノム編集のような編集効率の悪い場合に見落としが発生する。
  +  -> 2022-10-29 `threshold`を0.01, `OPTICS(min_samples=0.001)`として1%の点変異を検出できることを確認した

+ [x] Ayabe-Task1のbarcode31.fq(intact flow/wt)において、連続するCの部位に欠失（`C`CCCCC）がありmutatedとなっていた。
  + シークエンスエラーに起因するものと考えられるため、補正する。
+ [x] `screen_diffloci`：連続する塩基において、**deletionのみ**有意差を厳しくする。

## コンセンサス

+ [x] ~~consensusが同じクラスタにおいて、小さいクラスタを一番大きなクラスタのほうにマージする~~
+ [x] **コンセンサスのアレルにfloxやalbinoを含める** (DONE: 2022-10-09)
  + AyabeTask1のようにInsertionのなかにSubがある場合、controlアレルにマッピングしているとSubの検出ができないため
+ [x] コントロールにマッピングしたCSSPLITと、クラスタリングの結果を統合

## レポート

+ [x] (2022-11-01) report/BAM/{sample_name}のようにsample_nameごとに出力ディレクトリを分けたほうがわかりやすい
  + report/.igvjsにおいては{sample_name}/`{control_name}_control.bam(.bai)`という形式で保存した

+ [x] レポートのディレクトリ構造を作る
  + report/HTML
  + report/FASTA
  + report/BAM
  + report/VCF

+ [x] HTMLファイルのコンセンサス
+ [x] SAMファイルの配列からマイクロホモロジーを取り除く `report_bam.remove_microhomology`
+ [x] FASTAファイル
+ [x] (test) `report_af.summary_allele`
+ [x] (test) `report_af.plot`
  + [x] モジュールのimport文を最初に持ってくる
+ [x] (test) `report_af.all_allele`
+ [x] アリル割合のplot_alleles.pyの草稿作成（2022-03-11）-> `read_af`モジュールに移動
+ [x] 8以上のアレル数にも対応

## DAJIN singlemode

+ [x] 実装

## DAJIN view

+ [x] (2022-11-01)標準エラー出力を捨てる
+ [x] ブラウザが自動で起動するようにする
+ [x] controlを含める

+ [x] `igv.js`で各アレルの代表的なリード20本程度を可視化する
  + [x] (core/main.py) `igv.js`用に各アレルから20本のリードを抽出し、`report/.igvjs`に保存
  + [x] (batchmode/report.py) `index.html`を`DAJINResults/{name}/BAM/igvjs/`に保存して、`DAJINResults/{name}/BAM/.igvjs/`内にあるBAMとFASTA(reference用)を可視化する
    + [x] GENOMEがある場合
    + [x] GENOMEがない場合
  + [x] `DAJIN2 view -n/--name`で起動
+ [x] 20本ではなくて100本くらいにする？

## DAJIN gui
+ [x] (2022-11-01)標準エラー出力を捨てる
+ [x] ポート番号が毎回変わるようにする

## その他
