# Todo


## 2021-07-26
+ inversionアレルの変異情報が消えてしまう → 完成が先なので, 一旦無視します. 
+ paddingの”=”のあつかいはどうしよう. MatchかDeletionか？

---
## midsconvリファクタリング

+ [x] テストケースの記載
+ [x] テストケースの実行

### メモ
```sh
#* controlアレル以外のアレルに対する構造多型はすべて"controlの構造多型"としたいため, ミスマッチスコアを高くする
if (allele !~ /(control)|(wt)/) {
  c_of[1] = padD(gsub(/[MIDS]/, "", c_of[1]))
}
```
+ 上記のコードの必要性があるのか, もう一度検討する必要がありそうです.
  + 仮にref=control, que=inversionでマッピングしたら, すべてMatchという結果になってしまう.
    + その場合はflagで逆位になっているリードを`V`に変換したほうが良い
    + そのコードを実装した場合, 下記のコードはいらなくなる.
  + 仮にref=control, que=large_deletionでマッピングしたら, 
  + おそらく, Stx2 BC25のように大型欠失がある場合, Deletionと判定されてしまうのを妨げたいようです. ただし現在のコンセンサスはすべてcontrolを基準にするため, 必要がないと思います.

+ midsconvよりもsamTomidsvのほうが良い. 