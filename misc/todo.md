# ToDO

## 開発一般

+ [x] GitHub Actionsによるテスト自動化

## 前処理

+ [x] 大型欠失とInversionは断端のマイクロホモロジーがあるので、配列とQualityを元にして結合する -> `midsv v0.7.0`で実装
+ [ ] もとのSAMファイルの配列からマイクロホモロジーを取り除く

## マッピング
+ [x] mappyによるマッピング


## クラスタリング

## コンセンサス

+ [ ] コンセンサスの結果、同じクラスタリングの結果になったら統合する
+ [ ] VCFファイル
+ [ ] HTMLファイル
  + [ ] CRISPResso2のように割合もわかるとヘテロのときに嬉しい
+ [x] コントロールにマッピングしたCSSPLITと、クラスタリングの結果を統合
<!-- + [x] difflociのところだけ、コンセンサスコール。 -->
+ [ ] (test) call_percentage
+ [ ] (test) join_listdicts

## レポート
+ [x] FASTAファイル
+ [x] HTMLファイル
+ [ ] VCFファイル
+ [ ] BAMファイル
+ [x] (test) `report_af.summary_allele`
+ [x] (test) `report_af.plot`
+ [x] (test) `report_af.all_allele`

## 可視化
+ [x] アリル割合のplot_alleles.pyの草稿作成（2022-03-11）-> `read_af`モジュールに移動
+ [x] 8以上のアレル数にも対応
+ [ ] colorblind-friendly
+ [ ] `igv.js`で各アレルの代表的なリード20本程度を可視化する
  + [ ] `igv.js`は`npx http-server -a localhost -o . -s`を使うとローカルのBAMファイルを可視化できる

## その他

+ [ ] サンプルのみ（コントロールなし）のときにも（可能な限り）アレルを分離できるようにしたい
+ [ ] `CSSPLIT`は名が体を表し切れていないので改名したい
+ [ ] ターゲットアンプリコンではないとき（染色体にマッピングするとき）、MIDSV変換が巨大な文字列となってしまうときの対応
+ [ ] Batchモードの搭載（というかBatchを標準にする）

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
