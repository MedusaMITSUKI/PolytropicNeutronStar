# PolytropicNeutronStar
 ポリトロープモデルを用いた中性子星の状態方程式をTOV方程式で解く．

## 詳細
[YouTubeの動画](http://example.com)をご視聴ください．

## 必要なツール
- Fortranコンパイラ (gfortran等)
- Gnuplot 

## 使い方
1. Fortranのプログラムをコンパイルします．
   ```bash
   $ gfortran PolytropicNeutronStar.f90
   ```
   すると実行ファイル ``` a.out ``` が生成されます．

2. プログラムの実行
   ```a.out```を実行します．
   ```bash
   $ ./a.out
   ```
   すると計算ログファイル```polytrope.dat```とグラフ描画プログラム(Gnuplot)```polytrope.plt```が生成され，Gnuplotで計算結果が描画されます．

 ## パラメータについて
 Fortranプログラム内のGamma1の値を2.0 ~ 4.0の値に書き換えてみましょう．
 ポリトロピック指数(n, Gamma := 1 + 1/n)を変えることで，状態方程式の硬さが変化し，星の半径-質量のグラフの形が変わることが理解できます．
