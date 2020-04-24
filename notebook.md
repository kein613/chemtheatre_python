
# Part6 仮説検定

## Chap.0 全体の流れ
<p>Part 6では、仮説検定を行う。仮説検定とは、データから読み取れる仮定が、統計的に妥当かどうかを検証する操作のことである。<br>
具体的には、データを統計的に分析し、帰無仮説<i>H</i><sub>0</sub>が確率的に有り得そうか否かを判定する。この際、判断の基準として有意水準（&alpha;）を設定する。</p>
<p>一般的な検定の手順としては、まず、帰無仮説<i>H</i><sub>0</sub>を立て、有意水準（&alpha;）を設定する。次に、帰無仮説<i>H</i><sub>0</sub>の下で、データから求められる検定統計量を取る確率（<i>p</i>値）を求める。その後、その確率（<i>p</i>値）を有意水準（&alpha;）と比較する。<i>p</i>値が有意水準より小さい場合、確率的に帰無仮説<i>H</i><sub>0</sub>は滅多にありえないとして棄却し、対立仮説<i>H</i><sub>1</sub>を採択する。対して、。<i>p</i>値が有意水準より大きい場合、確率的に帰無仮説が正しいとして、帰無仮説<i>H</i><sub>0</sub>を採択する。</p>
<p>今回はChemTHEATREのいくつかのデータを用いて、仮説検定の概要を掴むことを目指す。</p>

<img src="../img/img07.SVG" alt="img07" style="zoom:80%;" />




## Chap.1 ライブラリの読み込み


```python
%matplotlib inline
import numpy as np
from scipy import stats
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
```

<p>まず最初にライブラリの読み込みを行う。今回利用するライブラリはいずれも、Anacondaにすべてインストールされているものである。</p>
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
    <td align="left">検定の計算に利用</td>
    <td align="left"><a href=https://www.scipy.org>https://www.scipy.org</a></td>
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
  <tr>
    <td align="left">seaborn</td>
    <td align="left">統計データ可視化ライブラリ</td>
    <td align="left">データの可視化に利用</td>
    <td align="left"><a href=https://seaborn.pydata.org>https://seaborn.pydata.org</a></td>
  </tr>
</table>



## Chap.2 データの読み込み

<p>ライブラリの読み込みができたら、次はデータの読み込みを行う。今回取り扱うデータは、スナメリ（<i>Neophocaena phocaenoides</i>）とスジイルカ（<i>Stenella coeruleoalba</i>）のデータである。これらをChemTHEATREのSample Searchからmeasureddataとsamplesを検索し、ダウンロードしてこのノートブックファイルに読み込む。</p>
```python
data_file1 = "measureddata_20190826061813.tsv"    #変数に入力する文字列を、各自のmeasureddataのtsvファイル名に変更する
data_file2 = "measureddata_20190826061826.tsv"    #変数に入力する文字列を、各自のmeasureddataのtsvファイル名に変更する

data1, data2 = pd.read_csv(data_file1, delimiter="\t"), pd.read_csv(data_file2, delimiter="\t")

data1 = data1.drop(["ProjectID", "ScientificName", "RegisterDate", "UpdateDate"], axis=1)    #後でsamplesと結合する際に重複する列の削除
data2 = data2.drop(["ProjectID", "ScientificName", "RegisterDate", "UpdateDate"], axis=1)    #後でsamplesと結合する際に重複する列の削除
```


```python
sample_file1 = "samples_20190826061810.tsv"    #変数に入力する文字列を、各自のsamplesのtsvファイル名に変更する
sample_file2 = "samples_20190826061824.tsv"    #変数に入力する文字列を、各自のsamplesのtsvファイル名に変更する

sample1, sample2 = pd.read_csv(sample_file1, delimiter="\t"), pd.read_csv(sample_file2, delimiter="\t")
```



## Chap.3 データの下処理

<p>それでは、読み込んだデータの下処理から行う。最初にsamplesとmeasureddataを結合する。</p>
```python
df1, df2 = pd.merge(data1, sample1, on="SampleID"), pd.merge(data2, sample2, on="SampleID")
```

<p>続いて、必要なデータの抽出を行う。まず、このデータフレームの単位を統一するために、単位が[ng/g lipid]のデータのみを抽出する。<br>
次に、抽出したデータからさらに4種類の化学物質（ΣPCBs, ΣDDTs, ΣPBDEs, ΣCHLs）のデータを抽出する。</p>


```python
data_lipid_1 = df1[(df1["Unit"] == "ng/g lipid")]
data_lipid_1 = data_lipid_1[(data_lipid_1["ChemicalName"] == "ΣPCBs") | (data_lipid_1["ChemicalName"] == "ΣDDTs") | 
                            (data_lipid_1["ChemicalName"] == "ΣPBDEs") |  (data_lipid_1["ChemicalName"] == "ΣCHLs")]

data_lipid_2 = df2[(df2["Unit"] == "ng/g lipid") & df2["ChemicalName"].str.startswith("Σ")]
data_lipid_2 = data_lipid_2[(data_lipid_2["ChemicalName"] == "ΣPCBs") | (data_lipid_2["ChemicalName"] == "ΣDDTs") | 
                            (data_lipid_2["ChemicalName"] == "ΣPBDEs")|  (data_lipid_2["ChemicalName"] == "ΣCHLs")]
```

<p>抽出したスナメリ・スジイルカのデータは結合して、空白の列を削除しておく。</p>
```python
data_lipid = pd.concat([data_lipid_1, data_lipid_2], sort = True)
data_lipid = data_lipid.dropna(how="all", axis=1)
```

<p>データがひとまとめになったので、ここからは必要なデータを切り出す。</p>
```python
data_pcb = data_lipid[data_lipid["ChemicalName"] == "ΣPCBs"]    #ΣPCBs
data_chl = data_lipid[data_lipid["ChemicalName"] == "ΣCHLs"]    #ΣCHLs
data_ddt = data_lipid[data_lipid["ChemicalName"] == "ΣDDTs"]    #ΣDDTs
data_pbde = data_lipid[data_lipid["ChemicalName"] == "ΣPBDEs"]    #ΣPBDEs
```


```python
data_1 = data_lipid[data_lipid["ScientificName"] == "Stenella coeruleoalba"]    #スジイルカのデータ
data_1 = data_1[(data_1["ChemicalName"] == "ΣPCBs") | (data_1["ChemicalName"] == "ΣCHLs")]    #スジイルカのデータからΣPCBsとΣCHLsを抽出
data_1.loc[:, "LogValue"] = data_1["MeasuredValue"].apply(np.log10)    #対数を取る
data_1
```



<table border="1" class="dataframe" style="font-size: 0.8rem">
  <thead>
    <tr style="text-align: left;">
      <th></th>
      <th>Age</th>
      <th>ChemicalID</th>
      <th>Chemical-Name</th>
      <th>Collection-Area</th>
      <th>Collection-Country</th>
      <th>Collection-Day</th>
      <th>Collection-LatitudeFrom</th>
      <th>Collection-LatitudeTo</th>
      <th>Collection-LongitudeFrom</th>
      <th>...</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>269</th>
      <td>11.5 y.o.</td>
      <td>CH0000096</td>
      <td>ΣPCBs</td>
      <td>Taiji, Wakayama</td>
      <td>Japan</td>
      <td>19.0</td>
      <td>33.544167</td>
      <td>33.623058</td>
      <td>135.885000</td>
      <td>...</td>
    </tr>
    <tr>
      <th>298</th>
      <td>11.5 y.o.</td>
      <td>CH0000152</td>
      <td>ΣCHLs</td>
      <td>Taiji, Wakayama</td>
      <td>Japan</td>
      <td>19.0</td>
      <td>33.544167</td>
      <td>33.623058</td>
      <td>135.885000</td>
      <td>...</td>
    </tr>
    <tr>
      <th>358</th>
      <td>17.5 y.o.</td>
      <td>CH0000096</td>
      <td>ΣPCBs</td>
      <td>Taiji, Wakayama</td>
      <td>Japan</td>
      <td>19.0</td>
      <td>33.544167</td>
      <td>33.623058</td>
      <td>135.885000</td>
      <td>...</td>
    </tr>
    <tr>
      <th>387</th>
      <td>17.5 y.o.</td>
      <td>CH0000152</td>
      <td>ΣCHLs</td>
      <td>Taiji, Wakayama</td>
      <td>Japan</td>
      <td>19.0</td>
      <td>33.544167</td>
      <td>33.623058</td>
      <td>135.885000</td>
      <td>...</td>
    </tr>
    <tr>
      <th>447</th>
      <td>19.5 y.o.</td>
      <td>CH0000096</td>
      <td>ΣPCBs</td>
      <td>Taiji, Wakayama</td>
      <td>Japan</td>
      <td>15.0</td>
      <td>33.544167</td>
      <td>33.623058</td>
      <td>135.885000</td>
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
      <th>1404</th>
      <td>NaN</td>
      <td>CH0000152</td>
      <td>ΣCHLs</td>
      <td>Gogo-shima, Ehime</td>
      <td>Japan</td>
      <td>1.0</td>
      <td>33.864628</td>
      <td>33.929025</td>
      <td>132.652733</td>
      <td>...</td>
    </tr>
    <tr>
      <th>1464</th>
      <td>NaN</td>
      <td>CH0000096</td>
      <td>ΣPCBs</td>
      <td>Gogo-shima, Ehime</td>
      <td>Japan</td>
      <td>1.0</td>
      <td>33.864628</td>
      <td>33.929025</td>
      <td>132.652733</td>
      <td>...</td>
    </tr>
    <tr>
      <th>1490</th>
      <td>NaN</td>
      <td>CH0000152</td>
      <td>ΣCHLs</td>
      <td>Gogo-shima, Ehime</td>
      <td>Japan</td>
      <td>1.0</td>
      <td>33.864628</td>
      <td>33.929025</td>
      <td>132.652733</td>
      <td>...</td>
    </tr>
    <tr>
      <th>1496</th>
      <td>NaN</td>
      <td>CH0000096</td>
      <td>ΣPCBs</td>
      <td>Gogo-shima, Ehime</td>
      <td>Japan</td>
      <td>1.0</td>
      <td>33.864628</td>
      <td>33.929025</td>
      <td>132.652733</td>
      <td>...</td>
    </tr>
    <tr>
      <th>1523</th>
      <td>NaN</td>
      <td>CH0000152</td>
      <td>ΣCHLs</td>
      <td>Gogo-shima, Ehime</td>
      <td>Japan</td>
      <td>1.0</td>
      <td>33.864628</td>
      <td>33.929025</td>
      <td>132.652733</td>
      <td>...</td>
    </tr>
  </tbody>
</table>
<p>42 rows × 39 columns</p>

<p>さらに、相関関係を見るために測定した ΣPCBs, ΣDDTs, ΣCHLs, ΣPBDEsのデータを同じスナメリの個体（SampleID）ごとにまとめる。</p>

<img src="../img/img13.SVG" alt="img13" style="zoom:80%;" />



```python
data_corr = pd.DataFrame(index = ["ScientificName", "ΣPCBs", "ΣDDTs", "ΣCHLs", "ΣPBDEs"])    #出力するDataFrame
for irow in range(len(sample1)):
    sample_id = sample1.at[irow, "SampleID"]
    if sample_id in data_lipid["SampleID"].values:
        rowdata = pd.Series(index = ["ScientificName", "ΣPCBs", "ΣDDTs", "ΣCHLs", "ΣPBDEs"])    #入力する行のフレーム
        rowdata["ScientificName"] = sample1.at[irow, "ScientificName"]    #SampleIDの検索
        for chem in ["ΣPCBs", "ΣDDTs", "ΣCHLs", "ΣPBDEs"]:    
            row = data_lipid_1[(data_lipid['SampleID'] == sample_id) & 
                               (data_lipid['ChemicalName'] == chem)].reset_index(drop = True)
            if row.empty:
                pass    #入力する行データがない場合、処理をしない
            else:
                rowdata[chem] = row.at[0, "MeasuredValue"]    #DataFrameに行データを追加
        data_corr = data_corr.append(rowdata, ignore_index=True)
    else:
        pass
data_corr = data_corr.dropna(how='any').reset_index(drop=True)
data_corr
```

    C:\Users\masah\Anaconda3\lib\site-packages\ipykernel_launcher.py:9: UserWarning: Boolean Series key will be reindexed to match DataFrame index.
      if __name__ == '__main__':



<table border="1" class="dataframe" style="font-size: 0.8rem">
  <thead>
    <tr style="text-align: left;">
      <th></th>
      <th>ScientificName</th>
      <th>ΣCHLs</th>
      <th>ΣDDTs</th>
      <th>ΣPBDEs</th>
      <th>ΣPCBs</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>Neophocaena phocaenoides</td>
      <td>770.0</td>
      <td>68000.0</td>
      <td>170.0</td>
      <td>5700.0</td>
    </tr>
    <tr>
      <th>1</th>
      <td>Neophocaena phocaenoides</td>
      <td>1200.0</td>
      <td>140000.0</td>
      <td>120.0</td>
      <td>6500.0</td>
    </tr>
    <tr>
      <th>2</th>
      <td>Neophocaena phocaenoides</td>
      <td>1000.0</td>
      <td>140000.0</td>
      <td>86.0</td>
      <td>5500.0</td>
    </tr>
    <tr>
      <th>3</th>
      <td>Neophocaena phocaenoides</td>
      <td>950.0</td>
      <td>130000.0</td>
      <td>100.0</td>
      <td>5800.0</td>
    </tr>
    <tr>
      <th>4</th>
      <td>Neophocaena phocaenoides</td>
      <td>1400.0</td>
      <td>280000.0</td>
      <td>140.0</td>
      <td>11000.0</td>
    </tr>
    <tr>
      <th>5</th>
      <td>Neophocaena phocaenoides</td>
      <td>1200.0</td>
      <td>130000.0</td>
      <td>91.0</td>
      <td>5000.0</td>
    </tr>
    <tr>
      <th>6</th>
      <td>Neophocaena phocaenoides</td>
      <td>1500.0</td>
      <td>220000.0</td>
      <td>84.0</td>
      <td>9200.0</td>
    </tr>
    <tr>
      <th>7</th>
      <td>Neophocaena phocaenoides</td>
      <td>340.0</td>
      <td>63000.0</td>
      <td>740.0</td>
      <td>4700.0</td>
    </tr>
    <tr>
      <th>8</th>
      <td>Neophocaena phocaenoides</td>
      <td>400.0</td>
      <td>51000.0</td>
      <td>780.0</td>
      <td>7200.0</td>
    </tr>
    <tr>
      <th>9</th>
      <td>Neophocaena phocaenoides</td>
      <td>140.0</td>
      <td>26000.0</td>
      <td>230.0</td>
      <td>1400.0</td>
    </tr>
    <tr>
      <th>10</th>
      <td>Neophocaena phocaenoides</td>
      <td>1400.0</td>
      <td>260000.0</td>
      <td>980.0</td>
      <td>28000.0</td>
    </tr>
    <tr>
      <th>11</th>
      <td>Neophocaena phocaenoides</td>
      <td>1900.0</td>
      <td>260000.0</td>
      <td>840.0</td>
      <td>22000.0</td>
    </tr>
    <tr>
      <th>12</th>
      <td>Neophocaena phocaenoides</td>
      <td>160.0</td>
      <td>38000.0</td>
      <td>320.0</td>
      <td>3000.0</td>
    </tr>
    <tr>
      <th>13</th>
      <td>Neophocaena phocaenoides</td>
      <td>520.0</td>
      <td>140000.0</td>
      <td>480.0</td>
      <td>10000.0</td>
    </tr>
    <tr>
      <th>14</th>
      <td>Neophocaena phocaenoides</td>
      <td>430.0</td>
      <td>110000.0</td>
      <td>470.0</td>
      <td>18000.0</td>
    </tr>
    <tr>
      <th>15</th>
      <td>Neophocaena phocaenoides</td>
      <td>150.0</td>
      <td>9900.0</td>
      <td>280.0</td>
      <td>1900.0</td>
    </tr>
    <tr>
      <th>16</th>
      <td>Neophocaena phocaenoides</td>
      <td>1400.0</td>
      <td>95000.0</td>
      <td>380.0</td>
      <td>14000.0</td>
    </tr>
    <tr>
      <th>17</th>
      <td>Neophocaena phocaenoides</td>
      <td>210.0</td>
      <td>43000.0</td>
      <td>220.0</td>
      <td>6700.0</td>
    </tr>
    <tr>
      <th>18</th>
      <td>Neophocaena phocaenoides</td>
      <td>440.0</td>
      <td>73000.0</td>
      <td>290.0</td>
      <td>9300.0</td>
    </tr>
    <tr>
      <th>19</th>
      <td>Neophocaena phocaenoides</td>
      <td>410.0</td>
      <td>61000.0</td>
      <td>290.0</td>
      <td>7800.0</td>
    </tr>
    <tr>
      <th>20</th>
      <td>Neophocaena phocaenoides</td>
      <td>320.0</td>
      <td>60000.0</td>
      <td>250.0</td>
      <td>2700.0</td>
    </tr>
    <tr>
      <th>21</th>
      <td>Neophocaena phocaenoides</td>
      <td>150.0</td>
      <td>25000.0</td>
      <td>89.0</td>
      <td>1400.0</td>
    </tr>
    <tr>
      <th>22</th>
      <td>Neophocaena phocaenoides</td>
      <td>160.0</td>
      <td>19000.0</td>
      <td>71.0</td>
      <td>2100.0</td>
    </tr>
  </tbody>
</table>



<p>最後に今後の計算で扱いやすいように、取り扱うデータをndarray形式に変換しておく。</p>
```python
pcb_4_1 = np.array(data_1[data_1["ChemicalName"] == "ΣPCBs"]["MeasuredValue"])    # スジイルカのΣPCBsデータ
chl_4_1 = np.array(data_1[data_1["ChemicalName"] == "ΣCHLs"]["MeasuredValue"])    # スジイルカのΣCHLsデータ

n_4_2 = np.array(data_pcb[data_pcb["ScientificName"] == "Neophocaena phocaenoides"]["MeasuredValue"])    # スナメリのΣPCBsデータ
s_4_2 = np.array(data_pcb[data_pcb["ScientificName"] == "Stenella coeruleoalba"]["MeasuredValue"])    # スジイルカのΣPCBsデータ

pcb_5_1 = np.array(data_corr["ΣPCBs"])    #スナメリのΣPCBsデータ
ddt_5_1 = np.array(data_corr["ΣDDTs"])    #スナメリのΣDDTsデータ

chl_5_2 = np.array(data_corr["ΣCHLs"])    #スナメリのΣCHLsデータ
pbde_5_2 = np.array(data_corr["ΣPBDEs"])    #スナメリのΣPBDEsデータ
```



## Chap.4 2標本t検定

<p>　2つの独立した母集団があるとして、それぞれの母集団から抽出した標本の平均に差があるかどうかを検定することを「2標本t検定」という。<br>
例えば、異なるクラスのテストの結果の比較や投薬前後での被験者の血圧の比較などに使われる。今回はChap.3で準備したChemTHEATREのデータを用いてこれを行う。</p>
<p>　2標本t検定については、それぞれの母集団の特性によって、行う操作が変わってくる。まず、検討すべきは「データに対応がある」か否かである。<br>
データに対応があるというのは、先の例で言う投薬前後での血圧の変化のように、データが同一の対象から抽出された「対」となる場合である。<br>
この場合は、2データ間の「差」を取ることができるので、「2群間の差が0であること」を検定する。<br>
　「データに対応がない」場合は、それぞれの母集団がどのような性質を持っているかによって、下図のように行う検定が異なってくる。<br>
また、これらの母集団の各性質についても、下図のようにその都度検定を行って確かめる必要がある。</p>
<img src="../img/img14.SVG" alt="img14" style="zoom:80%;" />


### Sec.4-1 スジイルカのΣPCBs・ΣCHLsデータの検定


```python
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
sns.stripplot(x="ChemicalName", y="LogValue", data=data_1, color='black', ax=ax)
ax.set_ylabel("ng/g lipid")
plt.suptitle("Stenella coeruleoalba")
plt.show()
```


![png](output_06_01.png)

<p>まず、スジイルカのΣPCBs・ΣCHLsデータの2つを比較してみる。この2つのデータについては、上の出力結果からΣPCBsのデータのほうが全体的に高い傾向にあると直感的にわかるが、それが統計的に妥当かどうかを検定を用いて検証してみる。</p>
<p>最初に行うのはシャピロ・ウィルク検定である。これは母集団が正規分布に従っているものかどうかの検定である。したがって今回の、帰無仮説<i>H</i><sub>0</sub>は「標本が正規分布からサンプリングされた」とするので、対立仮説<i>H</i><sub>1</sub>は「標本が正規分布からサンプリングされていない」である。<br>また、有意水準は5%（&alpha; = 0.05）として検定を行う。</p>
<p>PythonではScipyのstats.shapiro関数を利用すれば、シャピロ・ウィルク検定を行うことができる。この関数は返り値として検定統計量<i>W</i>と<i>p</i>値が得られる。</p>
```python
stats.shapiro(pcb_4_1), stats.shapiro(chl_4_1)
```




    ((0.9746639132499695, 0.8323823809623718),
     (0.9715853929519653, 0.7678512930870056))



<p>上の出力結果より、ΣPCBsのデータもΣCHLsのデータも共に、<i>p</i>値（右の値）が有意水準の0.05より大きい。したがって、帰無仮説<i>H</i><sub>0</sub>が採択されるので、これらのデータは正規分布からサンプリングされたといえる。</p>
<p>シャピロ・ウィルク検定から、データの正規性が確認できたので、次にF検定を行う。これは、これら2つのデータの分散が等しいことを検定するものである。なので、今回の帰無仮説<i>H</i><sub>0</sub>は「これら2つのデータの標準偏差は等しい」とし、対立仮説<i>H</i><sub>1</sub>は「2つのデータの標準偏差は等しくない」である。また、今回も有意水準5%（&alpha; = 0.05）で検定する。</p>
```python
stats.bartlett(pcb_4_1, chl_4_1)
```




    BartlettResult(statistic=18.10172409550588, pvalue=2.094117301987653e-05)



<p>上の出力結果から検定の結果を読み取ると、<i>p</i>値が有意水準の0.05より遥かに小さい。つまり、このような帰無仮説<i>H</i><sub>0</sub>は、まず滅多に起きず棄却され、対立仮説<i>H</i><sub>1</sub>が採択される。なので、スジイルカのΣPCBsとΣCHLsのデータは等分散性を仮定することができない。</p>
<p>ここまでの検定から、スジイルカのΣPCBs・ΣCHLsデータは正規性はあるものの、等分散ではないということがわかった。したがって、この2つのデータの比較は、ウェルチのt検定を行う。</p>
<p>PythonではScipyのstats.ttest_ind関数を使うことで、対応のない2標本のt検定を行うことができる。ちなみにequal_varパラメータをTrueにすると等分散と仮定したスチューデントのt検定を、Falseにすると等分散と仮定しないウェルチのt検定を行うことができる。</p>
```python
stats.ttest_ind(pcb_4_1, chl_4_1, equal_var=False)
```




    Ttest_indResult(statistic=9.046390311758614, pvalue=2.3876258572758015e-09)



<p>stats.ttest_ind関数の出力結果は、statisticがt検定の検定統計量<i>t</i>で、pvalueが<i>p</i>値である。<br>
今回は、有意水準（&alpha; = 0.05$）であり、<i>p</i>値がこれより小さいので、帰無仮説<i>H</i><sub>0</sub>は棄却され、対立仮説<i>H</i><sub>1</sub>の「ΣPCBsとΣCHLsのデータの平均値には差がある」が採択される。</p>
<p>検定結果を確認する意味でも、箱ひげ図を出力してみる。箱ひげ図の出力はPart1で行ったものと同じである。</p>
```python
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
sns.stripplot(x="ChemicalName", y="LogValue", data=data_1, color='black', ax=ax)
sns.boxplot(x="ChemicalName", y="LogValue", data=data_1, ax=ax)
ax.set_ylabel("ng/g lipid")
plt.suptitle("Stenella coeruleoalba")
plt.show()
```


![png](output_06_02.png)


### Sec.4-2 スナメリとスジイルカのΣPCBsデータの検定


```python
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
sns.stripplot(x="ScientificName", y="MeasuredValue", data=data_pcb, color='black', ax=ax)
sns.boxplot(x="ScientificName", y="MeasuredValue", data=data_pcb, ax=ax)
ax.set_ylabel("log(ng/g lipid)")
plt.suptitle("ΣPCBs")
plt.show()
```


![png](output_06_03.png)

<p>続いて、スナメリとスジイルカのΣPCBsのデータの比較を行ってみる。上の出力結果からも、スジイルカの方がスナメリより全体的に値が高い傾向にあることがわかる。今度も、統計的検定からそれを実証していく。</p>
<p>まず、2群間の統計検定をする前に、これらのデータ群にそれぞれ正規性があるかを確認するためにシャピロ・ウィルク検定を行う。<br>
帰無仮説<i>H</i><sub>0</sub>はSec.4-1同様に「標本は正規分布からサンプリングされた」であり、有意水準は5%（&alpha; = 0.05）とする。</p>


```python
stats.shapiro(n_4_2), stats.shapiro(s_4_2)    #n_4_2はスナメリ、s_4_2はスジイルカのデータ
```




    ((0.8355685472488403, 0.0011936090886592865),
     (0.9746639132499695, 0.8323823809623718))



<p>上の出力結果から、スナメリのΣPCBsデータの方は、<i>p</i>値は有意水準（&alpha; = 0.05）より小さく、正規性が認められなかった。<br>
よって、この2つのデータセットの比較にはt検定ができないので、代わりにマン・ホイットニーのU検定をする。<br>
マン・ホイットニーのU検定は、正規性の仮定できない独立した2群間の代表値の差がないことを検定する。<br>
つまり今回の場合、帰無仮説<i>H</i><sub>0</sub>は「スナメリのΣPCBsとスジイルカのΣPCBsのデータの代表値には差がない」とし、対立仮説<i>H</i><sub>1</sub>は「スナメリのΣPCBsとスジイルカのΣPCBsのデータの代表値には差がある」とする。また、有意水準は5% (&alpha; = 0.05)とする。</p>
<p>Pythonでマン・ホイットニーのU検定を行うには、Scipyのstats.mannwhitneyu関数を用いる。この関数の返り値は、検定統計量<i>U</i>と<i>p</i>値である。</p>
```python
stats.mannwhitneyu(n_4_2, s_4_2, alternative='two-sided')
```




    MannwhitneyuResult(statistic=41.5, pvalue=1.7554602710055643e-06)



<p>上の出力結果から、<i>p</i>値が有意水準（&alpha; = 0.05）より小さいので、帰無仮説<i>H</i><sub>0</sub>は棄却され、対立仮説<i>H</i><sub>1</sub>が採択される。<br>
つまり、スナメリとスジイルカのデータには差があるということがわかった。</p>



## Chap.5  無相関検定

<p>それでは、次に2変数の相関性を検定してみる。2変数の相関関係については、一般に散布図での可視化と相関係数の大きさで判断されことが多い。<br>
 このような相関性の有無の確認については、無相関検定を行う。無相関検定とは、帰無仮説<i>H</i><sub>0</sub>を「相関がない」として、それが起こる確率（<i>p</i>値）を有意水準（&alpha;）と比較する検定である。</p>

### Sec.5-1 スナメリのΣPCBsとΣDDTsのデータの相関の検定
<p>それでは早速、スナメリのΣPCBsとΣDDTsのデータの相関について考察してみる。<br>
まずは、matplotlibで散布図を描画してデータの概形を確認し、Numpyのcorrcoef関数を利用して相関係数を出してみる。</p>


```python
ax = plt.axes()
ax.scatter(x = pcb_5_1, y = ddt_5_1)
ax.set_xlabel("ΣPCBs [ng/g lipid]")
ax.set_ylabel("ΣDDTs [ng/g lipid]")
ax.set_title("Neophocaena phocaenoides")
plt.show()
```


![png](output_06_04.png)



```python
np.corrcoef(pcb_5_1, ddt_5_1)[0,1]
```




    0.7250244461730119



<p>上の散布図と相関係数から、正の相関があると考えられる。</p>
<p>次に、無相関検定を行って検証していく。まず、帰無仮説<i>H</i><sub>0</sub>は「スナメリのΣPCBsとΣDDTsのデータには相関がない」とする。また、有意水準は5%（&alpha; = 0.05）とする。</p>
<p>Pythonでは、Scipyのpearsonr関数を利用すれば、無相関検定を行うことができる。返り値は相関係数と<i>p</i>値である。</p>
```python
stats.pearsonr(pcb_5_1, ddt_5_1)
```




    (0.7250244461730119, 9.089030418519821e-05)



<p>上の出力結果より、<i>p</i>値が有意水準（&alpha; = 0.05）よりはるかに小さいので、帰無仮説<i>H</i><sub>0</sub>は棄却され、対立仮説の「無相関ではない」が採択される。<br>したがって、相関係数は有意だと考えられるので、スナメリのΣPCBsとΣDDTsのデータには正の相関があると言える。</p>
### Sec.5-2 スナメリのΣPBDEsとΣCHLsのデータの相関の検定
<p>最後に、スナメリのΣPBDEsとΣCHLsのデータの相関について検定してみる。今回も最初に散布図と相関係数を出力して考察する。</p>
```python
ax = plt.axes()
ax.scatter(x = pbde_5_2, y = chl_5_2)
ax.set_xlabel("ΣPBDEs [ng/g lipid]")
ax.set_ylabel("ΣCHLs [ng/g lipid]")
ax.set_title("Neophocaena phocaenoides")
plt.show()
```


![png](output_06_05.png)



```python
np.corrcoef(pbde_5_2, chl_5_2)[0,1]
```




    0.16324825566753187



<p>上の出力結果から、相関係数がゼロに近く、散布図からも相関性が明白に出ず、なんとなく無相関な感じである。このような状況での相関関係の考察において、無相関検定は有効である。</p>
<p>今回の無相関検定は、帰無仮説<i>H</i><sub>0</sub>は「ΣPBDEsとΣCHLsは無相関である」とし、有意水準は5%（&alpha; = 0.05）とする。</p>
```python
#stats.pearsonrで相関係数も見ることができるのにnp.corrcoeを単独で実行しておく理由がわからない
#また、私は正規性・等分散検定をする必要はない (自分の中の仮説に則って行えば良い) と考えていますが、
#実行するのであれば無相関検定においてなんの説明もなくピアソンの相関係数を使うのは良いのでしょうか。

stats.pearsonr(pbde_5_2, chl_5_2)
```




    (0.16324825566753182, 0.45671573485971206)



<p>出力結果は上のとおりである。$\rm{p}$値を確認すると、有意水準（&alpha; = 0.05）より大きいので、帰無仮説<i>H</i><sub>0</sub>が採択される。<br>
したがって、ΣPBDEsとΣCHLsとは相関関係がない事がわかった。</p>