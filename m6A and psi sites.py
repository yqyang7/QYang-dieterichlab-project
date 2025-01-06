>>> import pandas as pd
>>> m6A_df = pd.read_csv("/prj/Qingyang/m6A_cov10.bed", sep=r'\s+', low_memory=False)
>>> print(m6A_df.head())
  chrom  chromStart  chromEnd name  score strand  thickStart  thickEnd  itemRgb  coverage  frequency
0     1       14378     14379    a     10      -       14378     14379  255,0,0        10        0.0
1     1       14381     14382    a     17      -       14381     14382  255,0,0        17        0.0
2     1       14382     14383    a     14      -       14382     14383  255,0,0        14        0.0
3     1       14383     14384    a     12      -       14383     14384  255,0,0        12        0.0
4     1       14385     14386    a     11      -       14385     14386  255,0,0        11        0.0
>>> psi_df = pd.read_csv("/prj/Qingyang/psi_cov10.bed", sep=r'\s+', low_memory=False)
>>> print(psi_df.head())
  chrom  chromStart  chromEnd   name  score strand  thickStart  thickEnd  itemRgb  coverage  frequency
0     1       14384     14385  17802     11      -       14384     14385  255,0,0        11       9.09
1     1       14387     14388  17802     15      -       14387     14388  255,0,0        15       0.00
2     1       14390     14391  17802     10      -       14390     14391  255,0,0        10      30.00
3     1       14393     14394  17802     14      -       14393     14394  255,0,0        14       7.14
4     1       14395     14396  17802     12      -       14395     14396  255,0,0        12       0.00
>>> m6A_filtered = m6A_df[m6A_df['frequency'] >= 50.0]
>>> print(m6A_filtered.head())
    chrom  chromStart  chromEnd name  score strand  thickStart  thickEnd  itemRgb  coverage  frequency
13      1       14414     14415    a    149      -       14414     14415  255,0,0       149      65.10
39      1       14516     14517    a    162      -       14516     14517  255,0,0       162      78.40
67      1       14637     14638    a    159      -       14637     14638  255,0,0       159      81.76
88      1       14729     14730    a    107      -       14729     14730  255,0,0       107      59.81
121     1       14959     14960    a     14      -       14959     14960  255,0,0        14      64.29
>>> psi_filtered = psi_df[psi_df['frequency'] >= 50.0]
>>> print(psi_filtered.head())
    chrom  chromStart  chromEnd   name  score strand  thickStart  thickEnd  itemRgb  coverage  frequency
71      1       14780     14781  17802     67      -       14780     14781  255,0,0        67      64.18
72      1       14781     14782  17802    109      -       14781     14782  255,0,0       109      55.05
326     1       16056     16057  17802     15      -       16056     16057  255,0,0        15      66.67
461     1       16590     16591  17802     14      -       16590     16591  255,0,0        14      78.57
481     1       16761     16762  17802     66      -       16761     16762  255,0,0        66      53.03
>>> from sklearn.cluster import DBSCAN
>>> import numpy as np
>>> m6A_filtered = m6A_df[m6A_df['frequency'] >= 50.0].copy()
>>> psi_filtered = psi_df[psi_df['frequency'] >= 50.0].copy()
>>> m6A_filtered['center'] = (m6A_filtered['chromStart'] + m6A_filtered['chromEnd']) / 2
>>> psi_filtered['center'] = (psi_filtered['chromStart'] + psi_filtered['chromEnd']) / 2
>>> def perform_clustering(data):
...     clustering = DBSCAN(eps=1000, min_samples=2).fit(data[['center']])
...     return clustering
>>> m6A_clustering = perform_clustering(m6A_filtered)
>>> psi_clustering = perform_clustering(psi_filtered)
>>> m6A_filtered['cluster'] = m6A_clustering.labels_
>>> psi_filtered['cluster'] = psi_clustering.labels_
>>> print("m6A Clustering results:")
m6A Clustering results:
>>> print(m6A_filtered[['chrom', 'center', 'cluster']].head())
    chrom   center  cluster
13      1  14414.5        0
39      1  14516.5        0
67      1  14637.5        0
88      1  14729.5        0
121     1  14959.5        0
>>> print("psi Clustering results:")
psi Clustering results:
>>> print(psi_filtered[['chrom', 'center', 'cluster']].head())
    chrom   center  cluster
71      1  14780.5        0
72      1  14781.5        0
326     1  16056.5        0
461     1  16590.5        0
481     1  16761.5        0
>>> plt.figure(figsize=(12, 6))
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
NameError: name 'plt' is not defined
>>> import matplotlib.pyplot as plt
>>> from sklearn.cluster import DBSCAN
>>> plt.figure(figsize=(12, 6))
<Figure size 1200x600 with 0 Axes>
>>> plt.subplot(1, 2, 1)
<Axes: >
>>> plt.scatter(m6A_filtered['center'], np.zeros_like(m6A_filtered['center']), c=m6A_filtered['cluster'], cmap='viridis', marker='o')
<matplotlib.collections.PathCollection object at 0x7f92b8b49a50>
>>> plt.title('m6A Cluster Distribution')
Text(0.5, 1.0, 'm6A Cluster Distribution')
>>> plt.xlabel('Chromosome Position (center)')
Text(0.5, 0, 'Chromosome Position (center)')
>>> plt.ylabel('Density')
Text(0, 0.5, 'Density')
>>> plt.colorbar(label='Cluster ID')
<matplotlib.colorbar.Colorbar object at 0x7f92b8b86050>
>>> plt.subplot(1, 2, 2)
<Axes: >
>>> plt.scatter(psi_filtered['center'], np.zeros_like(psi_filtered['center']), c=psi_filtered['cluster'], cmap='viridis', marker='o')
<matplotlib.collections.PathCollection object at 0x7f92b8399190>
>>> plt.title('psi Cluster Distribution')
Text(0.5, 1.0, 'psi Cluster Distribution')
>>> plt.xlabel('Chromosome Position (center)')
Text(0.5, 0, 'Chromosome Position (center)')
>>> plt.ylabel('Density')
Text(0, 0.5, 'Density')
>>> plt.colorbar(label='Cluster ID')
<matplotlib.colorbar.Colorbar object at 0x7f92b83f4050>
>>> plt.tight_layout()
>>> plt.show()
>>> m6A_cluster_sizes = m6A_filtered['cluster'].value_counts()
>>> psi_cluster_sizes = psi_filtered['cluster'].value_counts()
>>> print("m6A size of each cluster:")
m6A size of each cluster:
>>> print(m6A_cluster_sizes)
cluster
-1        7999
 4586     1234
 8389      418
 6358      275
 11956     211
          ... 
 8981        2
 11981       2
 11980       2
 11976       2
 11971       2
Name: count, Length: 12021, dtype: int64
>>> print("psi size of each cluster:")
psi size of each cluster:
>>> print(psi_cluster_sizes)
cluster
-1        10913
 7183       437
 3300       399
 0          363
 9064       349
          ...  
 11747        2
 11750        2
 11752        2
 5325         2
 5346         2
Name: count, Length: 16124, dtype: int64
>>> merged = pd.merge(m6A_filtered, psi_filtered, on='center', suffixes=('_m6A', '_psi'))
>>> print("m6A and psi site intersection:")
m6A and psi site intersection:
>>> print(merged[['chrom_m6A', 'center', 'cluster_m6A', 'cluster_psi']].head())
  chrom_m6A      center  cluster_m6A  cluster_psi
0         1  26170858.5          304          401
1         1  39756326.5          479          672
2         1  42700803.5          520          722
3         1  54053166.5          643          948
4         1  93909398.5          795         2974
