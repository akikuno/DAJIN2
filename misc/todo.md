# ToDO

## 開発一般


+ [ ] 引数チェック関数 → メイン関数 → レポート関数の3部構成とする
  + メイン関数はbatchでも使えるようにする
  + 引数チェックとレポートはbatchとsingleで異なる（可能性がある）
  + -> (2022-10-07) batch.py, single.pyのなかに引数チェック関数とレポート関数を入れる方式にする


## 前処理

+ [ ] (test) `format_inputs`

## マッピング

## クラスタリング

+ [ ] (test) `clustering`

## コンセンサス

+ [ ] **コンセンサスのアレルにfloxやalbinoを含める**
  + AyabeTask1のようにInsertionのなかにSubがある場合、controlアレルにマッピングしているとSubの検出ができないため
+ [ ] (test) call_percentage
+ [ ] (test) join_listdicts

## レポート
+ [ ] HTMLファイルで、CRISPResso2のように割合もわかるとヘテロのときに嬉しい
  + [ ] colorblind-friendly
+ [ ] VCFファイル
+ [ ] BAMファイル
  + [x] Tyrでマイナス鎖からプラス鎖に変換
  + [ ] XXXでプラス鎖の確認
  + [ ] (test) `report_bam`
+ [ ] `igv.js`で各アレルの代表的なリード20本程度を可視化する
  + [ ] `igv.js`は`npx http-server -a localhost -o . -s`を使うとローカルのBAMファイルを可視化できる → Pythonで完結できそう…

## その他

+ [ ] サンプルのみ（コントロールなし）のときにも（可能な限り）アレルを分離できるようにしたい
+ [ ] `CSSPLIT`は名が体を表し切れていないので改名したい
+ [ ] ターゲットアンプリコンではないとき（染色体にマッピングするとき）、MIDSV変換が巨大な文字列となってしまうときの対応

---

# MEMO

## tyr

+ アルビノ点変異は829番目の塩基

## AyabeTask1

+ left-loxpは`control, SV false`にある。例: 880d5333-7a52-4dad-95fb-b778f06c1e5f
  + left-loxpがInsertion (at 1735) として扱われている。 left-loxpは46塩基長（<50）のためSVではない。
+ right-loxpは`flox, SV false`にある。例：6f536ce8-32d6-41d5-8b09-66bef425b9d8
  + left-loxpがDeletionとして扱われる。 left-loxpは46塩基長（<50）のためSVではない。


---
# DONE

## 開発一般

+ [x] (2022-10-07) `.tempdir-{name}`ではなく`.tempdir/{name}`とする (2022-10-08)
+ [x] batchモードの搭載
+ [x] GitHub Actionsによるテスト自動化

+ [x] `is_control`の実装
  + batchモードの時にcontrolが計算済みであれば処理を省略する
  + controlのファイルサイズが同じで、midsvのjsonlがあればコントロールは処理済みとする
  + メイン関数内のcontrolとsampleの処理を分ける

## 前処理
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

+ [x] `screen_diffloci`：連続する塩基において、**deletionのみ**有意差を厳しくする。

## コンセンサス

+ [x] コントロールにマッピングしたCSSPLITと、クラスタリングの結果を統合

## レポート

+ [x] HTMLファイルのコンセンサス
+ [x] SAMファイルの配列からマイクロホモロジーを取り除く `report_bam.remove_microhomology`
+ [x] FASTAファイル
+ [x] (test) `report_af.summary_allele`
+ [x] (test) `report_af.plot`
  + [x] モジュールのimport文を最初に持ってくる
+ [x] (test) `report_af.all_allele`
+ [x] アリル割合のplot_alleles.pyの草稿作成（2022-03-11）-> `read_af`モジュールに移動
+ [x] 8以上のアレル数にも対応

## その他
