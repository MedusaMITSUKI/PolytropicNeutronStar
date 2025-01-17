# PolytropicNeutronStar
 Polytropeモデルの中性子星の状態方程式から星の半径と質量を求める．

## 詳細
[YouTubeの動画](https://youtu.be/aq2uTYEiGvg)をご視聴ください．

## 必要なツール
- Fortranコンパイラ (gfortran等)
- Gnuplot 

## 使い方
Fortranのプログラムをコンパイルし，実行すると
- ```polytrope.dat``` ... 計算ログファイル
- ```polytrope.plt``` ... グラフ描画用のGnuplotファイル

が生成され，Gnuplotに計算結果が描画されます．
計算ログファイルには，左から順に星の質量 [太陽質量]・半径 [km]・中心密度 [kg/m^3] の計算結果が記録されます．

 ## パラメータについて
 Fortranプログラム内のGamma1の値を2.0 ~ 4.0の値に書き換えてみましょう．
 ```math
 where:
   Gamma := 1 + 1/n
   n ... polytropic指数
 ```
 ポリトロピック指数を変えることで，状態方程式の硬さが変化し，星の半径-質量のグラフの形が変わることが理解できます．
