# ToDO

## 開発一般

+ [x] GitHub Actionsによるテスト自動化

## 前処理

+ [x] 大型欠失とInversionは断端のマイクロホモロジーがあるので、Qualityを元にして結合する -> `midsv v0.7.0`で実装

### mapping.py
+ [x] mappyによるマッピング


## クラスタリング

## コンセンサス

+ [ ] コンセンサスの結果、同じクラスタリングの結果になったら統合する
+ [x] コントロールにマッピングしたCSSPLITと、クラスタリングの結果を統合
<!-- + [x] difflociのところだけ、コンセンサスコール。 -->
+ [x] FASTAファイル
+ [ ] VCFファイル
+ [ ] HTMLファイル
  + [ ] CRISPResso2のように割合もわかるとヘテロのときに嬉しい

## 可視化
+ [x] アリル割合のplot_alleles.pyの草稿作成（2022-03-11）
+ [ ] colorblind-friendly
+ [ ] 8以上のアレル数にも対応
+ [ ] `igv.js`で各アレルの代表的なリード20本程度を可視化する

## その他

+ [ ] `CSSPLIT`は名が体を表し切れていないので改名したい
+ [ ] ターゲットアンプリコンではないとき（染色体にマッピングするとき）、MIDSV変換が巨大な文字列となってしまうときの対応


---
# DONE

+ [x] 入力ファイルのフォーマット (format_input.py)
  + [x] sample, controlのテキスト化
  + [x] alleleでmulti-fastaをsingle fastaにする
  + [x] alleleでwtかcontrolがない場合にはエラー終了
  + [x] alleleで重複配列がある場合にはエラー終了

### samToMIDS.shをPythonで書き換え中… (midsvパッケージ)
+ [x] insertion
+ [x] substitution
+ [x] deletion
+ [x] padding
+ [x] mix condition
+ [x] Large deletion
+ [x] Large inversion
+ [x] CSSPLITしか使っておらず、無駄にMIDSVとQSCOREの計算時間がかかっている…→MIDSVパッケージの引数で計算するかしないかのオプションを付けたほうが良さそう（midsv=False, cssplit=True, qscore=False, など）
