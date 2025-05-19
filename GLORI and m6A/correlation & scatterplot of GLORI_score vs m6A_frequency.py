>>> import pandas as pd
>>> data = pd.read_csv("/prj/Qingyang/m6A_cov10.bed", sep='\t', low_memory=False)
>>> print(data.head(10))
  chrom  chromStart  chromEnd name  score strand  thickStart  thickEnd  itemRgb  coverage  frequency
0     1       14378     14379    a     10      -       14378     14379  255,0,0        10       0.00
1     1       14381     14382    a     17      -       14381     14382  255,0,0        17       0.00
2     1       14382     14383    a     14      -       14382     14383  255,0,0        14       0.00
3     1       14383     14384    a     12      -       14383     14384  255,0,0        12       0.00
4     1       14385     14386    a     11      -       14385     14386  255,0,0        11       0.00
5     1       14397     14398    a     12      -       14397     14398  255,0,0        12      25.00
6     1       14399     14400    a     23      -       14399     14400  255,0,0        23       0.00
7     1       14400     14401    a     21      -       14400     14401  255,0,0        21       4.76
8     1       14403     14404    a    142      -       14403     14404  255,0,0       142       1.41
9     1       14404     14405    a     98      -       14404     14405  255,0,0        98       2.04
>>> glo = pd.read_csv("/prj/Qingyang/GLORI.chrALL.tsv", sep=r'\s+', low_memory=False)
>>> print(glo.head(10))
  chrom  chromStart  chromEnd name  score strand ref5mer
0     1     1600376   1600377  m6A   31.1      +   GGACC
1     1     1600473   1600474  m6A   15.4      +   GGACT
2     1     1600621   1600622  m6A   25.7      +   AGACT
3     1     1601001   1601002  m6A   25.7      +   GGACT
4     1     1601143   1601144  m6A   18.3      +   GGACC
5     1     1601203   1601204  m6A   56.6      +   AGACT
6     1     1601293   1601294  m6A   23.2      +   GAATT
7     1     1601515   1601516  m6A   45.5      +   AGACC
8     1     1601840   1601841  m6A   50.6      +   GGACA
9     1     1602008   1602009  m6A   24.4      +   AGACT
>>> merged_data = pd.merge(data, glo, on=['chrom', 'chromStart', 'chromEnd', 'strand'],how='inner',suffixes=('_data', '_glo'))
>>> result = merged_data[['chrom', 'chromStart', 'chromEnd', 'strand', 'score_glo', 'frequency']]
>>> result.columns = ['chrom', 'chromStart', 'chromEnd', 'strand', 'GLORI_score', 'm6A_frequency']
>>> result.to_csv("file_merged_data.bed", index=False, sep='\t')
>>> print(result)
       chrom  chromStart   chromEnd strand  GLORI_score  m6A_frequency
0          1      185037     185038      -         69.3          61.11
1          1      185232     185233      -         34.4          12.09
2          1      826938     826939      -         65.8          57.14
3          1      841506     841507      +         73.7          92.31
4          1      841558     841559      +         74.5          91.67
...      ...         ...        ...    ...          ...            ...
130673     X   155524439  155524440      -         24.5          18.75
130674     X   155524446  155524447      -         64.6          50.00
130675     X   155881357  155881358      +         18.2          21.43
130676     X   155942014  155942015      +         17.8          28.65
130677     X   155942139  155942140      +         13.2          19.61

[130678 rows x 6 columns]
>>> correlation = result['GLORI_score'].corr(result['m6A_frequency'])
>>> print(f"\correlation between GLORI_score and m6A_frequency：{correlation:.4f}")
\correlation between GLORI_score and m6A_frequency：0.9069
>>> import matplotlib
>>> matplotlib.use('Agg')
>>> import seaborn as sns
>>> import matplotlib.pyplot as plt
>>> sns.set(style="white")
>>> plt.figure(figsize=(8, 6))
<Figure size 800x600 with 0 Axes>
>>> correlation = result[['GLORI_score', 'm6A_frequency']].corr()
>>> print(correlation)
               GLORI_score  m6A_frequency
GLORI_score       1.000000       0.906932
m6A_frequency     0.906932       1.000000
>>> sns.heatmap(correlation, annot=True, cmap="coolwarm", fmt=".2f", linewidths=0.5)
<Axes: >
>>> plt.title("Correlation between GLORI_score and m6A_frequency")
Text(0.5, 1.0, 'Correlation between GLORI_score and m6A_frequency')
>>> plt.savefig('heatmap_new.png')
>>> plt.figure(figsize=(8, 6))
<Figure size 800x600 with 0 Axes>
>>> plt.figure(figsize=(8, 6))
<Figure size 800x600 with 0 Axes>
>>> sns.scatterplot(x='GLORI_score',y='m6A_frequency',data=result,color="blue",alpha=0.7,edgecolor=None)
<Axes: xlabel='GLORI_score', ylabel='m6A_frequency'>
>>> plt.title("Scatterplot of GLORI_score vs m6A_frequency", fontsize=14)
Text(0.5, 1.0, 'Scatterplot of GLORI_score vs m6A_frequency')
>>> plt.xlabel("GLORI_score", fontsize=12)
Text(0.5, 0, 'GLORI_score')
>>> plt.ylabel("m6A_frequency", fontsize=12)
Text(0, 0.5, 'm6A_frequency')
>>> sns.regplot(x='GLORI_score',y='m6A_frequency',data=result,scatter=False,color="red",line_kws={'linewidth': 2})
<Axes: title={'center': 'Scatterplot of GLORI_score vs m6A_frequency'}, xlabel='GLORI_score', ylabel='m6A_frequency'>
>>> plt.tight_layout()
>>> plt.show()
>>> plt.savefig('scatterplot of GLORI_score vs m6A_frequency.png')
>>> exit()
