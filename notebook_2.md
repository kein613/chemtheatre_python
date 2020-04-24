
# Part5 地図上へのプロット（後編）
<p>さて、前編の続きである。後編では、前編で出力したグラフを時系列ごとに並べて、パラパラ漫画の要領でアニメーションを作成する。</p>
<figure id="process">
<img src="../img/img (6).SVG" alt="全体図"width=75% height=75%>
</figure>



## Chap.6 ライブラリの読み込み


```python
import os
import glob
import itertools
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')     # matplotlibのグラフ保存の設定
import matplotlib.pyplot as plt
from cartopy import crs as ccrs
from PIL import Image, ImageDraw
```

<p>最初に利用するライブラリを読み込む事から始める。<br>今回は、JupyterNotebook上にグラフ表示するを必要がないので「%matplotlib inline」のマジックコマンドは入力しない。<br>
  その代わり、出力結果を保存する<sup><a href=#sup1>1</a></sup>ために、matplotlibの設定を調整する必要がある（7行目・コメント参照）。</p>
<table style="text-align: left; font-size: 0.8rem">
  <tr align="left">
    <th>ライブラリ</th>
    <th>概要</th>
    <th>今回の使用目的</th>
    <th>公式URL</th>
  </tr>
   <tr align="left">
    <td>os</td>
    <td>標準ライブラリ</td>
    <td>ディレクトリの操作</td>
    <td><a href=https://docs.python.org/ja/3/library/os.html>https://docs.python.org/ja/3/library/os.html</a></td>
  </tr>
  <tr align="left">
    <td>glob</td>
    <td>ファイル・ディレクトリへの<br>アクセスの標準ライブラリ</td>
    <td>ファイル一覧の取得</td>
    <td><a href=https://docs.python.org/ja/3/library/glob.html>https://docs.python.org/ja/3/library/glob.html</a></td>
  </tr>
  <tr align="left">
    <td>itertools</td>
    <td>イテレータを生成する<br>標準ライブラリ</td>
    <td>ループ処理を効率化するのに利用</td>
    <td><a href=https://docs.python.org/ja/3/library/itertools.html>https://docs.python.org/ja/3/library/itertools.html</a></td>
  </tr>
  <tr align="left">
    <td>NumPy</td>
    <td>数値計算ライブラリ</td>
    <td>統計処理上の数値計算に利用</td>
    <td><a href=https://www.numpy.org>https://www.numpy.org</a></td>
  </tr>
  <tr align="left">
    <td>pandas</td>
    <td>データ分析ライブラリ</td>
    <td>データ処理や整形に利用</td>
    <td><a href=https://pandas.pydata.org>https://pandas.pydata.org</a></td>
  </tr>
  <tr align="left">
    <td>Matplotlib</td>
    <td>グラフ描画ライブラリ</td>
    <td>データの可視化に利用</td>
    <td><a href=https://matplotlib.org>https://matplotlib.org</a></td>
  </tr>
  <tr align="left">
    <td>cartopy</td>
    <td>地図描画ライブラリ</td>
    <td>地図データの可視化に利用</td>
    <td><a href=https://scitools.org.uk/cartopy/docs/latest>https://scitools.org.uk/cartopy/docs/latest</a></td>
  </tr>
  <tr align="left">
    <td>Pillow<sup><a href=#sup2>2</a></sup></td>
    <td>画像処理ライブラリ</td>
    <td>GIFアニメーションの作成に利用</td>
    <td><a href=https://pillow.readthedocs.io/en/stable/>https://pillow.readthedocs.io/en/stable/</a></td>
  </tr>
</table>



## Chap.7 データの準備 

<p>アニメーションを作成する前に、データの準備を行う。今回取り扱うデータは、Part 5の前編のChap.3で作成・出力し、グラフを表示したものと同じデータである。</p>
```python
data_lipid = pd.read_csv("data.csv")
data_lipid
```



<table border="1" class="dataframe" style="font-size: 0.8rem">
  <thead>
    <tr style="text-align: left;">
      <th></th>
      <th>Unnamed: 0</th>
      <th>MeasuredID</th>
      <th>SampleID</th>
      <th>ChemicalID</th>
      <th>Chemical-Name</th>
      <th>...</th>
      <th>Measured-Value</th>
      <th>Alternative-Data</th>
      <th>ProjectID</th>
      <th>...</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>4064</td>
      <td>23738</td>
      <td>SAA002109</td>
      <td>CH0000096</td>
      <td>ΣPCBs</td>
      <td>...</td>
      <td>30000.0</td>
      <td>NaN</td>
      <td>PRA000036</td>
      <td>...</td>
    </tr>
    <tr>
      <th>1</th>
      <td>4065</td>
      <td>23792</td>
      <td>SAA002109</td>
      <td>CH0000033</td>
      <td>ΣDDTs</td>
      <td>...</td>
      <td>68000.0</td>
      <td>NaN</td>
      <td>PRA000036</td>
      <td>...</td>
    </tr>
    <tr>
      <th>2</th>
      <td>4069</td>
      <td>24008</td>
      <td>SAA002109</td>
      <td>CH0000152</td>
      <td>ΣCHLs</td>
      <td>...</td>
      <td>2700.0</td>
      <td>NaN</td>
      <td>PRA000036</td>
      <td>...</td>
    </tr>
    <tr>
      <th>3</th>
      <td>4075</td>
      <td>24382</td>
      <td>SAA002109</td>
      <td>CH0000146</td>
      <td>ΣHCHs</td>
      <td>...</td>
      <td>300.0</td>
      <td>NaN</td>
      <td>PRA000036</td>
      <td>...</td>
    </tr>
    <tr>
      <th>4</th>
      <td>4081</td>
      <td>23739</td>
      <td>SAA002110</td>
      <td>CH0000096</td>
      <td>ΣPCBs</td>
      <td>...</td>
      <td>30000.0</td>
      <td>NaN</td>
      <td>PRA000036</td>
      <td>...</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>475</th>
      <td>6211</td>
      <td>25508</td>
      <td>SAA001987</td>
      <td>CH0000146</td>
      <td>ΣHCHs</td>
      <td>...</td>
      <td>380.0</td>
      <td>NaN</td>
      <td>PRA000032</td>
      <td>...</td>
    </tr>
    <tr>
      <th>476</th>
      <td>6277</td>
      <td>25425</td>
      <td>SAA001988</td>
      <td>CH0000096</td>
      <td>ΣPCBs</td>
      <td>...</td>
      <td>36000.0</td>
      <td>NaN</td>
      <td>PRA000032</td>
      <td>...</td>
    </tr>
    <tr>
      <th>477</th>
      <td>6278</td>
      <td>25446</td>
      <td>SAA001988</td>
      <td>CH0000033</td>
      <td>ΣDDTs</td>
      <td>...</td>
      <td>67000.0</td>
      <td>NaN</td>
      <td>PRA000032</td>
      <td>...</td>
    </tr>
    <tr>
      <th>478</th>
      <td>6279</td>
      <td>25467</td>
      <td>SAA001988</td>
      <td>CH0000152</td>
      <td>ΣCHLs</td>
      <td>...</td>
      <td>11000.0</td>
      <td>NaN</td>
      <td>PRA000032</td>
      <td>...</td>
    </tr>
    <tr>
      <th>479</th>
      <td>6281</td>
      <td>25509</td>
      <td>SAA001988</td>
      <td>CH0000146</td>
      <td>ΣHCHs</td>
      <td>...</td>
      <td>520.0</td>
      <td>NaN</td>
      <td>PRA000032</td>
      <td>...</td>
    </tr>
  </tbody>
</table>
<p>480 rows × 40 columns</p>


<p>ここで、アニメーションの出力に利用するので、いつからいつまでのデータが有るかを確認する。最小値と最大値はminメソッドとmaxメソッドを利用する。</p>
```python
start, end = data_lipid["CollectionYear"].min(), data_lipid["CollectionYear"].max()
start, end
```




    (1978, 2006)



<p>また、どのような生物・化学物質のデータを含むのかも確認しておく。これはuniqueメソッドを利用する。</p>
```python
lipid_species = data_lipid["ScientificName"].unique()
lipid_chemicals = data_lipid["ChemicalName"].unique()
lipid_species, lipid_chemicals
```




    (array(['Peponocephala electra', 'Neophocaena phocaenoides',
            'Sousa chinensis', 'Stenella coeruleoalba'], dtype=object),
     array(['ΣPCBs', 'ΣDDTs', 'ΣCHLs', 'ΣHCHs'], dtype=object))



## Chap.8 アニメーション化する

### Sec.8-1 時系列ごとにプロットする
<p>Chap.4・Chap.5で出力した画像は、すべての時間のデータを重ね合わせて出力していた。このSec.ではそのデータを時系列ごとにバラして出力する。</p>
<p>まず、保存先の準備をしておく。</p>
```python
save_dir = "fig"
os.mkdir(save_dir)   #フォルダの作成
```


    ---------------------------------------------------------------------------
    
    FileExistsError                           Traceback (most recent call last)
    
    <ipython-input-17-d07c45681d29> in <module>
          1 save_dir = "fig"
    ----> 2 os.mkdir(save_dir)   #フォルダの作成


    FileExistsError: [WinError 183] 既に存在するファイルを作成することはできません。: 'fig'


<p>保存先の準備ができたら、グラフの出力をする。</p>
```python
#このコード以下のコードについてについてなんの説明もなし？理解できるとは思えません。

for year in range(start, end+1):
    df_list = []
    for k1, k2 in itertools.product(lipid_chemicals, lipid_species):
        df_list.append(data_lipid[(data_lipid["CollectionYear"] == year) & (data_lipid["ChemicalName"] == k1) & (data_lipid["ScientificName"] == k2)])
    fig = plt.figure(figsize=(16, 9))
    ax = [0]*4
    rate = [10,100,10,10]
    for i in range(4): 
        ax[i] = fig.add_subplot(2, 2, i+1, projection=ccrs.PlateCarree())
        ax[i].coastlines()
        ax[i].scatter(x = np.array(df_list[4*i+0]["CollectionLongitudeFrom"]), y = np.array(df_list[4*i+0]["CollectionLatitudeFrom"]),
                      s=np.array(df_list[4*i+0]["MeasuredValue"])/rate[i], c = "red", alpha=0.3)
        ax[i].scatter(x = np.array(df_list[4*i+1]["CollectionLongitudeFrom"]), y = np.array(df_list[4*i+1]["CollectionLatitudeFrom"]),
                      s=np.array(df_list[4*i+1]["MeasuredValue"])/rate[i], c = "blue", alpha=0.3)
        ax[i].scatter(x = np.array(df_list[4*i+2]["CollectionLongitudeFrom"]), y = np.array(df_list[4*i+2]["CollectionLatitudeFrom"]),
                      s=np.array(df_list[4*i+2]["MeasuredValue"])/rate[i], c = "yellow", alpha=0.3)
        ax[i].scatter(x = np.array(df_list[4*i+3]["CollectionLongitudeFrom"]), y = np.array(df_list[4*i+3]["CollectionLatitudeFrom"]),
                      s=np.array(df_list[4*i+3]["MeasuredValue"])/rate[i], c = "green", alpha=0.3)
        ax[i].set_xlim(90,180)
        ax[i].set_ylim(15, 60)
        plt.title(lipid_chemicals[i])
    plt.suptitle(year)
    plt.savefig(os.path.join(save_dir, str(year)+".png"))
    plt.close()
```

### Sec.8-2 統合してアニメーション出力する
<p>Sec.8-1で、年ごとのグラフを出力したので、あとはパラパラ漫画の要領で画像を重ねて一つのGIFファイル<sup><a href=#sup5>5</a></sup>にする。</p>
<p>まず、出力した画像ファイルの一覧をリストで取得する。</p>
```python
f_list = glob.glob(os.path.join(save_dir, "*"))
```

<p>取得したファイルを順に読み込み、GIFアニメーションとしてまとめて出力する。画像の読み込み・書き出しはPillowを利用する。<br>
なお、アニメーションの速さは、saveメソッドのdurationというパラメーターでミリ秒単位での調整ができる。</p>


```python
images = []
for f in f_list:
    img = Image.open(f)
    images.append(img)
images[0].save(os.path.join(save_dir, "animation.gif"), save_all=True, append_images=images[1:], optimize=False, duration=500, loop=0)
```



## 脚注

<p><sup id=sup1>1</sup>Jupyterは本来サーバで動作するアプリケーションであり、JupyterNotebookは自分のPCに立ち上げたサーバにアクセスしている。<br>Jupyterは外部からの接続が前提のサーバアプリのため、基本的に出力結果をディスプレイに表示して、サーバには保存しないように設計されている。<br>そのため画像を出力保存する際は、表示するディスプレイを設定する代わりにAGG（Anti-Grain Geometry）という画像出力ライブラリを利用するように設定する必要がある。</p>
<p><sup id=sup2>2</sup>Pillowは、PILというPythonの画像ライブラリから派生したライブラリである。ソースコードを受け継いでいるのでコード内ではPILと表記されている。</p>
