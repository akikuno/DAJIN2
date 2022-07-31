# ToDO

## 開発一般

+ [x] GitHub Actionsによるテスト自動化

## 前処理

### mapping.py
+ [x] mappyによるマッピング

### samToMIDS.shをPythonで書き換え中… (midsvパッケージ)
+ [x] insertion
+ [x] substitution
+ [x] deletion
+ [x] padding
+ [x] mix condition
+ [x] Large deletion
+ [x] Large inversion

+ **CSSPLITしか使っておらず、無駄にMIDSVとQSCOREの計算時間がかかっている…→MIDSVパッケージの引数で計算するかしないかのオプションを付けたほうが良さそう（midsv=False, cssplit=True, qscore=False, など）**

## クラスタリング


## アリル割合の可視化
+ [x] plot_alleles.pyの草稿作成（2022-03-11）
+ [ ] colorblind-friendly
+ [ ] 8以上のアレル数にも対応


---
# DONE

+ [x] 入力ファイルのフォーマット (format_input.py)
  + [x] sample, controlのテキスト化
  + [x] alleleでmulti-fastaをsingle fastaにする
  + [x] alleleでwtかcontrolがない場合にはエラー終了
  + [x] alleleで重複配列がある場合にはエラー終了
