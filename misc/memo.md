# memo

# tyr

+ アルビノ点変異は829番目の塩基です

# controlのあつかい
+ controlはもしかしてwtアレルとcontrolアレルでしかつかわないかも？

## ライブラリとスクリプトの違い
+ テスト対象であるのがライブラリ
  + ライブラリの実態はシェル関数
+ テストしたライブラリを呼び出して使用するのがスクリプト

cat .DAJIN_temp/midsmask/barcode31_control.csv  | awk -F, '{print $(200+1)}' | sort | uniq -c

cat .DAJIN_temp/midsmask/barcode31_control.csv  | awk -F, '{print $(828+1)}' | sort | uniq -c
cat .DAJIN_temp/mids/barcode31_control.csv  | awk -F, '{print $(829+1)}' | sort | uniq -c
cat .DAJIN_temp/midsmask/barcode31_control.csv  | awk -F, '{print $(830+1)}' | sort | uniq -c

cat .DAJIN_temp/midsmask/barcode31_albino.csv  | awk -F, '{print $(828+1)}' | sort | uniq -c
cat .DAJIN_temp/midsmask/barcode31_albino.csv  | awk -F, '{print $(829+1)}' | sort | uniq -c
cat .DAJIN_temp/midsmask/barcode31_albino.csv  | awk -F, '{print $(830+1)}' | sort | uniq -c
