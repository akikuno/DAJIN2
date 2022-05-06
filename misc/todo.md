# ToDO

## 開発一般

+ [ ] GitHub Workflowによるテスト自動化
## 前処理

+ mapping.py
  + [ ] mappyによるマッピング

+ samToMIDS.shをPythonで書き換え中… (test_mids.py)
  + [x] insertion
  + [x] substitution
  + [x] deletion
  + [ ] padding
  + [ ] mix condition
  + [ ] Large deletion
  + [ ] Large inversion

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
