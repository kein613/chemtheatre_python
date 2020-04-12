
# Part 2 統計的推定



## Chap.0 全体の流れ

<p>Part2では、統計的推定を行う。自然科学では、調査対象（自然）は直接全て調査することができないので、標本調査が一般的である。<br>
その際、得られた標本（母集団の一部）から、未知の母集団全体の分布を推定するのが、統計的推定である。</p>
<figure id="process">
<img src="../img/img (2).SVG" alt="全体図" width=75% height=75%>
</figure>



## Chap.1 ライブラリの読み込み


```python
%matplotlib inline
import numpy as np
from scipy import stats
import math
import pandas as pd
from matplotlib import pyplot as plt
```



最初に、必要なライブラリの読み込みから始める。１行目は、Part 1同様にmatplotlibをJupyter Notebook内で表示するためのマジックコマンドである。

２行目以降が、今回利用するライブラリである。これらのライブラリのうち、mathライブラリはPythonに標準で組み込まれている。また、それ以外のライブラリはAnacondaにインストールされている。

ちなみに、mathとNumpyは機能が類似しているが、mathは標準のPythonに組み込まれている関数で、Numpyは複雑な数値計算を効率的にする拡張モジュールであり、その役割は異なる。



<table style="font-size: 0.8rem">
  <tr>
    <th>ライブラリ</th>
    <th>概要</th>
    <th>今回の使用目的</th>
    <th>公式URL</th>
  </tr>
  <tr>
    <td align="left">NumPy</td>
    <td align="left">数値計算ライブラリ</td>
    <td align="left">統計処理上の数値計算に利用</td>
    <td align="left"><a href=https://www.numpy.org>https://www.numpy.org</a></td>
  </tr>
  <tr>
    <td align="left">Scipy</td>
    <td align="left">科学計算ライブラリ</td>
    <td align="left">統計的推定の計算に利用</td>
    <td align="left"><a href=https://www.scipy.org>https://www.scipy.org</a></td>
  </tr>
  <tr>
    <td align="left">math</td>
    <td align="left">標準の数値計算ライブラリ</td>
    <td align="left">平方根などのかんたんな計算に利用</td>
    <td align="left"><a href=https://docs.python.org/ja/3/library/math.html>https://docs.python.org/ja/3/library/math.html</a></td>
  </tr>
  <tr>
    <td align="left">pandas</td>
    <td align="left">データ分析ライブラリ</td>
    <td align="left">データ読み込みや整形に利用</td>
    <td align="left"><a href=https://pandas.pydata.org>https://pandas.pydata.org</a></td>
  </tr>
  <tr>
    <td align="left">Matplotlib</td>
    <td align="left">グラフ描画ライブラリ</td>
    <td align="left">データの可視化に利用</td>
    <td align="left"><a href=https://matplotlib.org>https://matplotlib.org</a></td>
  </tr>
</table>



## Chap.2 データの読み込み

<p>今回はカツオ（<i>Katsuwonus pelamis</i>）のデータを利用する。Part 0のChap.4を参照し、ChemTHEATREのSample Searchから、カツオのサンプルデータを計測データをダウンロードする。<br>
ダウンロードできたら、このノートブックファイルのあるフォルダにmeasureddataとsamplesのデータを移動する。その後Anacondaを起動し直した後に、Part 1同様にpandasのread_csv関数を利用して、計測データと試料データの双方を読み込む。</p>




```python
data_file = "measureddata_20190930045953.tsv"    #変数に入力する文字列を、各自のmeasureddataのtsvファイル名に変更する
chem = pd.read_csv(data_file, delimiter="\t")
chem = chem.drop(["ProjectID", "ScientificName", "RegisterDate", "UpdateDate"], axis=1)    #後でsamplesと結合する際に重複する列の削除

sample_file = "samples_20190930045950.tsv"    #変数に入力する文字列を、各自のsamplesのtsvファイル名に変更する
sample = pd.read_csv(sample_file, delimiter="\t")
```

<p>pythonのようにプログラムでファイルを読み込んだり、加工したりした際は、想定したとおりにファイルが読み込めているか確認する癖付けをしておいたほうが良い。ちなみにJupyter Notebookの場合、変数名のみ入力すると、その変数の中身がOutに表示されるので便利である。</p>


```python
chem
```



<table border="1" class="dataframe" style="font-size: 0.6rem">
  <thead>
    <tr style="text-align: left;">
      <th></th>
      <th>MeasuredID</th>
      <th>SampleID</th>
      <th>ChemicalID</th>
      <th>ChemicalName</th>
      <th>ExperimentID</th>
      <th>MeasuredValue</th>
      <th>AlternativeData</th>
      <th>Unit</th>
      <th>Remarks</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>1</td>
      <td>SAA000001</td>
      <td>CH0000096</td>
      <td>ΣPCBs</td>
      <td>EXA000001</td>
      <td>6.659795</td>
      <td>NaN</td>
      <td>ng/g wet</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>1</th>
      <td>2</td>
      <td>SAA000002</td>
      <td>CH0000096</td>
      <td>ΣPCBs</td>
      <td>EXA000001</td>
      <td>9.778107</td>
      <td>NaN</td>
      <td>ng/g wet</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>2</th>
      <td>3</td>
      <td>SAA000003</td>
      <td>CH0000096</td>
      <td>ΣPCBs</td>
      <td>EXA000001</td>
      <td>5.494933</td>
      <td>NaN</td>
      <td>ng/g wet</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>3</th>
      <td>4</td>
      <td>SAA000004</td>
      <td>CH0000096</td>
      <td>ΣPCBs</td>
      <td>EXA000001</td>
      <td>7.354636</td>
      <td>NaN</td>
      <td>ng/g wet</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>4</th>
      <td>5</td>
      <td>SAA000005</td>
      <td>CH0000096</td>
      <td>ΣPCBs</td>
      <td>EXA000001</td>
      <td>9.390950</td>
      <td>NaN</td>
      <td>ng/g wet</td>
      <td>NaN</td>
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
    </tr>
    <tr>
      <th>74</th>
      <td>75</td>
      <td>SAA000082</td>
      <td>CH0000096</td>
      <td>ΣPCBs</td>
      <td>EXA000001</td>
      <td>3.321208</td>
      <td>NaN</td>
      <td>ng/g wet</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>75</th>
      <td>76</td>
      <td>SAA000083</td>
      <td>CH0000096</td>
      <td>ΣPCBs</td>
      <td>EXA000001</td>
      <td>3.285111</td>
      <td>NaN</td>
      <td>ng/g wet</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>76</th>
      <td>77</td>
      <td>SAA000084</td>
      <td>CH0000096</td>
      <td>ΣPCBs</td>
      <td>EXA000001</td>
      <td>0.454249</td>
      <td>NaN</td>
      <td>ng/g wet</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>77</th>
      <td>78</td>
      <td>SAA000085</td>
      <td>CH0000096</td>
      <td>ΣPCBs</td>
      <td>EXA000001</td>
      <td>0.100000</td>
      <td>&lt;1.00E-1</td>
      <td>ng/g wet</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>78</th>
      <td>79</td>
      <td>SAA000086</td>
      <td>CH0000096</td>
      <td>ΣPCBs</td>
      <td>EXA000001</td>
      <td>0.702224</td>
      <td>NaN</td>
      <td>ng/g wet</td>
      <td>NaN</td>
    </tr>
  </tbody>
</table>
<p>79 rows × 9 columns</p>



```python
sample
```

<table border="1" class="dataframe" style="font-size: 0.6rem">
  <thead>
    <tr style="text-align: left;">
      <th></th>
      <th>ProjectID</th>
      <th>SampleID</th>
      <th>SampleType</th>
      <th>TaxonomyID</th>
      <th>UniqCodeType</th>
      <th>UniqCode</th>
      <th>SampleName</th>
      <th>ScientificName</th>
      <th>CommonName</th>
      <th>CollectionYear</th>
      <th>...</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>PRA000001</td>
      <td>SAA000001</td>
      <td>ST008</td>
      <td>8226</td>
      <td>es-BANK</td>
      <td>EF00564</td>
      <td>NaN</td>
      <td>Katsuwonus pelamis</td>
      <td>Skipjack tuna</td>
      <td>1998</td>
      <td>...</td>
    </tr>
    <tr>
      <th>1</th>
      <td>PRA000001</td>
      <td>SAA000002</td>
      <td>ST008</td>
      <td>8226</td>
      <td>es-BANK</td>
      <td>EF00565</td>
      <td>NaN</td>
      <td>Katsuwonus pelamis</td>
      <td>Skipjack tuna</td>
      <td>1998</td>
      <td>...</td>
    </tr>
    <tr>
      <th>2</th>
      <td>PRA000001</td>
      <td>SAA000003</td>
      <td>ST008</td>
      <td>8226</td>
      <td>es-BANK</td>
      <td>EF00566</td>
      <td>NaN</td>
      <td>Katsuwonus pelamis</td>
      <td>Skipjack tuna</td>
      <td>1998</td>
      <td>...</td>
    </tr>
    <tr>
      <th>3</th>
      <td>PRA000001</td>
      <td>SAA000004</td>
      <td>ST008</td>
      <td>8226</td>
      <td>es-BANK</td>
      <td>EF00567</td>
      <td>NaN</td>
      <td>Katsuwonus pelamis</td>
      <td>Skipjack tuna</td>
      <td>1998</td>
      <td>...</td>
    </tr>
    <tr>
      <th>4</th>
      <td>PRA000001</td>
      <td>SAA000005</td>
      <td>ST008</td>
      <td>8226</td>
      <td>es-BANK</td>
      <td>EF00568</td>
      <td>NaN</td>
      <td>Katsuwonus pelamis</td>
      <td>Skipjack tuna</td>
      <td>1998</td>
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
      <td>...</td>
    </tr>
    <tr>
      <th>74</th>
      <td>PRA000001</td>
      <td>SAA000082</td>
      <td>ST008</td>
      <td>8226</td>
      <td>es-BANK</td>
      <td>EF00616</td>
      <td>NaN</td>
      <td>Katsuwonus pelamis</td>
      <td>Skipjack tuna</td>
      <td>1999</td>
      <td>...</td>
    </tr>
    <tr>
      <th>75</th>
      <td>PRA000001</td>
      <td>SAA000083</td>
      <td>ST008</td>
      <td>8226</td>
      <td>es-BANK</td>
      <td>EF00617</td>
      <td>NaN</td>
      <td>Katsuwonus pelamis</td>
      <td>Skipjack tuna</td>
      <td>1999</td>
      <td>...</td>
    </tr>
    <tr>
      <th>76</th>
      <td>PRA000001</td>
      <td>SAA000084</td>
      <td>ST008</td>
      <td>8226</td>
      <td>es-BANK</td>
      <td>EF00619</td>
      <td>NaN</td>
      <td>Katsuwonus pelamis</td>
      <td>Skipjack tuna</td>
      <td>1999</td>
      <td>...</td>
    </tr>
    <tr>
      <th>77</th>
      <td>PRA000001</td>
      <td>SAA000085</td>
      <td>ST008</td>
      <td>8226</td>
      <td>es-BANK</td>
      <td>EF00620</td>
      <td>NaN</td>
      <td>Katsuwonus pelamis</td>
      <td>Skipjack tuna</td>
      <td>1999</td>
      <td>...</td>
    </tr>
    <tr>
      <th>78</th>
      <td>PRA000001</td>
      <td>SAA000086</td>
      <td>ST008</td>
      <td>8226</td>
      <td>es-BANK</td>
      <td>EF00621</td>
      <td>NaN</td>
      <td>Katsuwonus pelamis</td>
      <td>Skipjack tuna</td>
      <td>1999</td>
      <td>...</td>
    </tr>
  </tbody>
</table>
<p>79 rows × 66 columns</p>



## Chap.3 データの下処理

<p>データの読み込みが完了したら、次はデータの下処理を行う。<br>
   まず、２つに分かれているデータ（chemとsample）を統合し、必要なデータのみ抽出する。今回は、カツオのΣPCBのデータを利用したいので、"ChemicalName"列の値が"ΣPCB"のデータのみを抽出する。</p>




```python
df = pd.merge(chem, sample, on="SampleID")
data = df[df["ChemicalName"] == "ΣPCBs"]
```



<p>続いて、計測データの単位が異なっているかどうかを確認する。Part 1のようにデータの単位が異なっていると単純に比較や統合ができないからである。</p>

```python
data["Unit"].unique()
```




    array(['ng/g wet'], dtype=object)

​	

<p>pandasのuniqueメソッドを利用すると、そのデータフレーム内に含まれる値の一覧を見ることができる。ここで、"Unit"列に含まれる値の一覧を出力してみると、"ng/g wet"のみである事がわかるので、今回は、単位によるデータの分割は不要である。</p>
<p>最後は、N/Aしかない列を削除して、データの下処理は完了である。</p>

```python
data = data.dropna(how='all', axis=1)
data
```



<table border="1" class="dataframe" style="font-size: 0.6rem">
  <thead>
    <tr style="text-align: left;">
      <th></th>
      <th>MeasuredID</th>
      <th>SampleID</th>
      <th>ChemicalID</th>
      <th>ChemicalName</th>
      <th>ExperimentID</th>
      <th>MeasuredValue</th>
      <th>AlternativeData</th>
      <th>Unit</th>
      <th>ProjectID</th>
      <th>SampleType</th>
      <th>...</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>1</td>
      <td>SAA000001</td>
      <td>CH0000096</td>
      <td>ΣPCBs</td>
      <td>EXA000001</td>
      <td>6.659795</td>
      <td>NaN</td>
      <td>ng/g wet</td>
      <td>PRA000001</td>
      <td>ST008</td>
      <td>...</td>
    </tr>
    <tr>
      <th>1</th>
      <td>2</td>
      <td>SAA000002</td>
      <td>CH0000096</td>
      <td>ΣPCBs</td>
      <td>EXA000001</td>
      <td>9.778107</td>
      <td>NaN</td>
      <td>ng/g wet</td>
      <td>PRA000001</td>
      <td>ST008</td>
      <td>...</td>
    </tr>
    <tr>
      <th>2</th>
      <td>3</td>
      <td>SAA000003</td>
      <td>CH0000096</td>
      <td>ΣPCBs</td>
      <td>EXA000001</td>
      <td>5.494933</td>
      <td>NaN</td>
      <td>ng/g wet</td>
      <td>PRA000001</td>
      <td>ST008</td>
      <td>...</td>
    </tr>
    <tr>
      <th>3</th>
      <td>4</td>
      <td>SAA000004</td>
      <td>CH0000096</td>
      <td>ΣPCBs</td>
      <td>EXA000001</td>
      <td>7.354636</td>
      <td>NaN</td>
      <td>ng/g wet</td>
      <td>PRA000001</td>
      <td>ST008</td>
      <td>...</td>
    </tr>
    <tr>
      <th>4</th>
      <td>5</td>
      <td>SAA000005</td>
      <td>CH0000096</td>
      <td>ΣPCBs</td>
      <td>EXA000001</td>
      <td>9.390950</td>
      <td>NaN</td>
      <td>ng/g wet</td>
      <td>PRA000001</td>
      <td>ST008</td>
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
      <td>...</td>
    </tr>
    <tr>
      <th>74</th>
      <td>75</td>
      <td>SAA000082</td>
      <td>CH0000096</td>
      <td>ΣPCBs</td>
      <td>EXA000001</td>
      <td>3.321208</td>
      <td>NaN</td>
      <td>ng/g wet</td>
      <td>PRA000001</td>
      <td>ST008</td>
      <td>...</td>
    </tr>
    <tr>
      <th>75</th>
      <td>76</td>
      <td>SAA000083</td>
      <td>CH0000096</td>
      <td>ΣPCBs</td>
      <td>EXA000001</td>
      <td>3.285111</td>
      <td>NaN</td>
      <td>ng/g wet</td>
      <td>PRA000001</td>
      <td>ST008</td>
      <td>...</td>
    </tr>
    <tr>
      <th>76</th>
      <td>77</td>
      <td>SAA000084</td>
      <td>CH0000096</td>
      <td>ΣPCBs</td>
      <td>EXA000001</td>
      <td>0.454249</td>
      <td>NaN</td>
      <td>ng/g wet</td>
      <td>PRA000001</td>
      <td>ST008</td>
      <td>...</td>
    </tr>
    <tr>
      <th>77</th>
      <td>78</td>
      <td>SAA000085</td>
      <td>CH0000096</td>
      <td>ΣPCBs</td>
      <td>EXA000001</td>
      <td>0.100000</td>
      <td>&lt;1.00E-1</td>
      <td>ng/g wet</td>
      <td>PRA000001</td>
      <td>ST008</td>
      <td>...</td>
    </tr>
    <tr>
      <th>78</th>
      <td>79</td>
      <td>SAA000086</td>
      <td>CH0000096</td>
      <td>ΣPCBs</td>
      <td>EXA000001</td>
      <td>0.702224</td>
      <td>NaN</td>
      <td>ng/g wet</td>
      <td>PRA000001</td>
      <td>ST008</td>
      <td>...</td>
    </tr>
  </tbody>
</table>
<p>79 rows × 35 columns</p>

## Chap.4 点推定

<p>標本から母集団を推測する、統計的推定のうち、ピンポイントで値を推定するのが点推定である。ここでは、カツオから検出されたΣPCB濃度の標本から、母集団（採集地域での個体全体）のΣPCB濃度を推定してみる。</p>
<p>まず、年ごとに計算しその変化の推移を見るために、何年のデータが含まれているかを確認する。</p>

```python
data['CollectionYear'].unique()
```




    array([1998, 1997, 1999, 2001], dtype=int64)



<p>上のuniqueメソッドから1997～1999年の3年間と2001年のデータがデータセットに含まれていることがわかった。ここでまずは、1997年のデータを取り出してみる。この際、今後の計算が楽になるように取り出したデータをNumpyのndarray<sup><a href=#sup1>1</a></sup>形式に変更しておく。</p>


```python
pcb_1997 = np.array(data[data['CollectionYear']==1997]["MeasuredValue"]) # 1997年の測定値のみを抽出
pcb_1997
```




    array([ 10.72603788,   9.22208078,   7.59790835,  30.95079465,
            15.27462553,  14.15719633,  13.28955903,  14.87712806,
             9.86650189,  18.26554514,   3.39951845,   6.58172781,
            12.43564814,   6.1948639 ,   6.41605666,   4.98827291,
            12.36669815,  31.17955551,   8.16184346,   4.60893266,
            36.85826409,  52.99841724,  39.22500351,  53.92302702,
            69.4308048 ,  73.97686479, 125.3887794 ,  45.39974771,
            54.12726127,  39.77794045, 101.2736126 ,  38.06220403,
           126.8301693 ,  70.25308435,  31.24246301,  21.3958656 ,
            41.85726522,  30.91112132,  81.12597135,  10.76755148,
            24.20442213,  24.57497594,  14.84353549,  59.53687389,
            52.78443082,   8.4644697 ,   4.15293758,   3.31957452,
             4.51832675,   6.98373973])



<p>同様に1998年と1999年のデータも抽出する。</p>

```python
pcb_1998 = np.array(data[data['CollectionYear']==1998]["MeasuredValue"]) # 1998年の測定値のみを抽出
pcb_1999 = np.array(data[data['CollectionYear']==1999]["MeasuredValue"]) # 1999年の測定値のみを抽出
```



<p>ここで、平均と分散の不偏推定量を算出する。<br>
まず、平均の不偏推定量（$\hat{\mu}$）だが、これは標本平均（$\overline{X}$）の期待値が母平均と等しいことを利用する。（下式参照）<br>
$$ E \left(\overline X \right) = E\left(\frac{1}{n} \sum_{i=1}^{n} \left(x_i\right)\right) = \frac{1}{n}\sum_{i=1}^{n} E\left(x_i\right) = \frac{1}{n} \times n\mu = \mu \\ \therefore \hat\mu = \overline{x} $$</p>




```python
s_mean_1997 = np.mean(pcb_1997)
s_mean_1997
```




    31.775384007760003



<p>同様に、分散の不偏推定量を算出する。<br>
このとき標本分散（$S^2$）の期待値は、母分散（$\sigma^2$）と同じ値は取らず、代わりに不偏分散（$s^2$）を求める必要があることに注意する。<br>
$$\hat\sigma^2 \neq S^2 = \frac{1}{n} \sum_{i=1}^{n} \left( x_i - \overline X \right) \\ \hat\sigma^2 = s^2 = \frac{1}{n-1} \sum_{i=1}^{n} \left( x_i - \overline X \right)$$</p>
<p>なお、Numpyのvar関数は、どちらの分散も算出することができ、不偏分散はddof=1のパラメータで出力される。ただし、デフォルトではddof=0の標本分散が出力させるので注意が必要である。<br>
$$\mathrm{np.var}\left(x_1 \ldots x_n, \mathrm{ddof=0}\right): S^2 = \frac{1}{n} \sum_{i=1}^{n} \left( x_i - \overline X \right) \\
\mathrm{np.var}\left(x_1 \ldots x_n, \mathrm{ddof=1}\right): \hat\sigma^2 = \frac{1}{n-1} \sum_{i=1}^{n} \left( x_i - \overline X \right)$$</p>




```python
u_var_1997 = np.var(pcb_1997, ddof=1)
u_var_1997 
```




    942.8421749786518



<p>同様に、1998年と1999年の平均と分散の不偏推定量を算出する。</p>

```python
s_mean_1998, s_mean_1999 = np.mean(pcb_1998), np.mean(pcb_1999)
u_var_1998, u_var_1999 = np.var(pcb_1998, ddof=1), np.var(pcb_1999, ddof=1)
```



<p>ここで、求めた代表値について整理する。まず、</p>

```python
s_mean_1997, s_mean_1998, s_mean_1999
```




    (31.775384007760003, 17.493267312533337, 30.583242522000003)




```python
u_var_1997, u_var_1998, u_var_1999
```




    (942.8421749786518, 240.2211176248311, 1386.7753819003349)



## Chap.5 区間推定と信頼区間

<p>Chap.5では、Chap.4で求めた点推定とは異なり、母平均や母分散を統計的に一定の範囲で推定する区間推定をする。</p>

### Sec.5-1 母平均の区間推定



<p>まず、区間推定をする前に、各年のデータセットのデータ数を調べる。データ数は、pythonに標準で実装されているlen関数を利用すれば、算出できる。</p>

```python
n_1997 = len(pcb_1997)
n_1997
```




    50




```python
n_1998, n_1999 = len(pcb_1998), len(pcb_1999)
n_1998, n_1999
```




    (15, 13)



<p>上記から、1997年～1999年の各年のデータセットのデータ総数がわかった。</p>
<p>この内、1997年のデータセットは、$n = 50$と大標本であり、1998年・1999年のデータセットは、それぞれ$n = \left\{ \begin{array}{ll}15 & \left( 1998 \right) \\ 13 & \left( 1999 \right) \end{array} \right.$で、小標本である。<br>したがって、このあとの区間推定の処理が少々異なることに注意する必要がある。</p>
<p>まず1997年のデータセットから、母平均の区間推定をする。この場合、母分散未知で大標本($n > 30$)なので、中心極限定理から標本平均($\overline X $)は正規分布 $N\left(  \mu ,  \frac{s^2}{n} \right)$を近似することができる。なので、母平均を信頼度（$\alpha$）で区間推定すると、信頼区間は下式のようになる。</p>
$\overline X - z_\frac{\alpha}{2} \sqrt{\frac{s^2}{n}} < \mu < \overline X - z_\frac{\alpha}{2} \sqrt{\frac{s^2}{n}} $
<p>pythonではScipyのstas.norm.interval()を利用すると、平均(loc)・標準偏差(scale)の正規分布でalpha×100%となる範囲を、中央値を中心として取得できる。<br>
ここで、信頼度（$\alpha = 0.95$）で信頼区間を算出する。</p>




```python
m_interval_1997 = stats.norm.interval(alpha=0.95, loc=s_mean_1997, scale=math.sqrt(pcb_1997.var(ddof=1)/n_1997))
m_interval_1997
```




    (23.26434483549182, 40.28642318002819)



<p>次に、1998年、1999年のデータセットについて母平均の区間推定をする。これらは、母分散が未知で、小標本($n\leq 30$)である。この場合、母平均$\mu$は、正規分布$ N \left( \mu , \frac{s^2}{n} \right)$ではなく、自由度（$n-1$）のt分布を利用する。なので、母平均を信頼度（$\alpha$）で区間推定すると、信頼区間は下式のようになる。</p>
$\overline X - t_\frac{\alpha}{2}\left(n-1\right)\sqrt{\frac{s^2}{n}} < \mu < \overline X + t_\frac{\alpha}{2}\left(n-1\right)\sqrt{\frac{s^2}{n}} $
<p>pythonでは、stats.t.interval()を利用すると、平均(loc)・標準偏差(scale)・自由度(df)のt分布でalpha×100%となる範囲を、中央値を中心として取得できる。<br>
ここでは、信頼度（$\alpha = 0.95 $）で信頼区間を算出する。</p>




```python
m_interval_1998 = stats.t.interval(alpha=0.95, df=n_1998-1, loc=s_mean_1998, scale=math.sqrt(pcb_1998.var(ddof=1)/n_1998))
m_interval_1999 = stats.t.interval(alpha=0.95, df=n_1999-1, loc=s_mean_1999, scale=math.sqrt(pcb_1999.var(ddof=1)/n_1999))
```




```python
m_interval_1997, m_interval_1998, m_interval_1999
```




    ((23.26434483549182, 40.28642318002819),
     (8.910169386248537, 26.076365238818138),
     (8.079678286109523, 53.086806757890486))



<p>なお、95%信頼区間とは、母平均が95%の確率でその範囲にあるということを表している。つまり、信頼度（$\alpha$）を小さくすると、母平均が信頼区間に含まれる確率が小さくなると同時に、信頼区間は狭くなる。</p>

```python
stats.norm.interval(alpha=0.9, loc=s_mean_1997, scale=math.sqrt(pcb_1997.var(ddof=1)/n_1997))
```




    (24.63269477364296, 38.91807324187704)





### Sec.5-2 母分散の区間推定



<p>次に、母分散の区間推定をする。母分散（$\sigma^2$）の区間推定では、$\frac{\left(n-1\right)s^2}{\sigma^2}$が、自由度$(n-1)$の$\chi^2$分布に従うことを利用する。</p>
$\chi_\frac{\alpha}{2}\left(n-1\right) \leq \frac{\left( n-1 \right)s^2}{\sigma^2} \leq \chi_{1-\frac{\alpha}{2}}\left(n-1\right)$
<p>まず、自由度（$n-1$）の$\chi^2$分布のパーセント点（$\chi_\frac{\alpha}{2}\left(n-1\right), \chi_{1-\frac{\alpha}{2}}\left(n-1\right)$）を求める。ここでは、信頼度0.95で計算する。<br>
なお、pythonではScipyのstats.chi2.interval()で、自由度（df）のalpha×100%となる範囲が取得できる。</p>




```python
chi_025_1997, chi_975_1997 = stats.chi2.interval(alpha=0.95, df=n_1997-1)
chi_025_1997, chi_975_1997
```




    (31.554916462667137, 70.22241356643451)



<p>続いて、信頼区間を求める。導出には、以下の式を参考にする。</p>
$\frac{\left(n-1\right)s^2}{\chi_\frac{\alpha}{2}\left(n-1\right)} \leq \sigma^2 \leq \frac{\left(n-1\right)s^2}{\chi_{1-\frac{\alpha}{2}}\left(n-1\right)}$


```python
v_interval_1997 = (n_1997 - 1)*np.var(pcb_1997, ddof=1) / chi_975_1997, (n_1997 - 1)*np.var(pcb_1997, ddof=1) / chi_025_1997
v_interval_1997
```




    (657.8991553778869, 1464.0909168183869)



<p>同様に、1998年、1999年のデータセットに関しても、分散の区間推定をする。なお、母平均の区間推定とは異なり、$\frac{\left(n-1\right)s^2}{\sigma^2}$の分布は、標本サイズに関わらず、$\chi^2$分布に従う。</p>

```python
chi_025_1998, chi_975_1998 = stats.chi2.interval(alpha=0.95, df=n_1998-1)
chi_025_1999, chi_975_1999 = stats.chi2.interval(alpha=0.95, df=n_1999-1)
```




```python
v_interval_1998 = (n_1998 - 1)*np.var(pcb_1998, ddof=1) / chi_975_1998, (n_1998 - 1)*np.var(pcb_1998, ddof=1) / chi_025_1998
v_interval_1999 = (n_1999 - 1)*np.var(pcb_1999, ddof=1) / chi_975_1999, (n_1999 - 1)*np.var(pcb_1999, ddof=1) / chi_025_1999
```




```python
v_interval_1997, v_interval_1998, v_interval_1999
```




    ((657.8991553778869, 1464.0909168183869),
     (128.76076176378118, 597.4878836139195),
     (713.0969734866349, 3778.8609867211235))



```python
chi_025_1997, chi_975_1997 = stats.chi2.interval(alpha=0.9, df=n_1997-1)
(n_1997 - 1)*np.var(pcb_1997, ddof=1) / chi_975_1997, (n_1997 - 1)*np.var(pcb_1997, ddof=1) / chi_025_1997
```




    (696.4155490924242, 1361.5929987004467)



## Chap.6 推定結果の可視化

<p>それでは、Chap.4・Chap.5で推定した母平均をグラフに可視化する。</p>
<p>まず、Chap.4で点推定した母平均の値を時系列にまとめる。</p>

```python
x_list = [1997, 1998, 1999]
y_list = [s_mean_1997, s_mean_1998, s_mean_1999]
```



<p>次に、Chap.5で推定した、信頼度95%の母平均の信頼区間も時系列にまとめる。</p>

```python
interval_list = []
interval_list.append(m_interval_1997)
interval_list.append(m_interval_1998)
interval_list.append(m_interval_1999)
interval_list
```




    [(23.26434483549182, 40.28642318002819),
     (8.910169386248537, 26.076365238818138),
     (8.079678286109523, 53.086806757890486)]



<p>母平均の95%信頼区間は、このままでは可視化に利用できないので、信頼区間の幅を求める。</p>

```python
interval_list = np.array(interval_list).T[1] - y_list
x_list, y_list, interval_list
```




    ([1997, 1998, 1999],
     [31.775384007760003, 17.493267312533337, 30.583242522000003],
     array([ 8.51103917,  8.58309793, 22.50356424]))



<p>最後に、matplotlibで可視化する。信頼区間は、エラーバーで表示する。matplotlibでエラーバーを表示する際は、errorbarメソッドを利用する。<br>
このメソッドでは、X軸の値（ここでは年）、Y軸の値（ここでは点推定の母平均）、エラーバーの長さ（ここでは信頼区間の幅）を指定する。</p>




```python
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.errorbar(x=x_list, y=y_list, yerr=interval_list, fmt='o-', capsize=4, ecolor='red')
plt.xticks(x_list)
ax.set_title("Katsuwonus pelamis")
ax.set_ylabel("ΣPCBs [ng/g wet]")
plt.show()
```

![png](output_59_0.png)




## 脚注
<p><sup id=sup1>1</sup>Numpyでのn次行列を格納するデータ形式。</p>
