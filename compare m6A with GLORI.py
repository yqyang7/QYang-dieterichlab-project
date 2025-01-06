>>> import matplotlib
>>> matplotlib.use('Agg')
>>> import seaborn as sns
>>> import matplotlib.pyplot as plt
>>> import pandas as pd
>>> data = pd.read_csv("/prj/Qingyang/m6A_cov10.bed", sep=r'\s+', low_memory=False)
>>> print(data.head())
  chrom  chromStart  chromEnd name  score strand  thickStart  thickEnd  itemRgb  coverage  frequency
0     1       14378     14379    a     10      -       14378     14379  255,0,0        10        0.0
1     1       14381     14382    a     17      -       14381     14382  255,0,0        17        0.0
2     1       14382     14383    a     14      -       14382     14383  255,0,0        14        0.0
3     1       14383     14384    a     12      -       14383     14384  255,0,0        12        0.0
4     1       14385     14386    a     11      -       14385     14386  255,0,0        11        0.0
>>> glo = pd.read_csv("/prj/Qingyang/GLORI.chrALL.tsv", sep=r'\s+', low_memory=False)
>>> print(glo.head())
  chrom  chromStart  chromEnd name  score strand ref5mer
0     1     1600376   1600377  m6A   31.1      +   GGACC
1     1     1600473   1600474  m6A   15.4      +   GGACT
2     1     1600621   1600622  m6A   25.7      +   AGACT
3     1     1601001   1601002  m6A   25.7      +   GGACT
4     1     1601143   1601144  m6A   18.3      +   GGACC
>>> columns_to_compare = ["chrom","chromStart","chromEnd","strand"]
>>> data_subset = data[columns_to_compare]
>>> glo_subset = glo[columns_to_compare]
>>> sametrans = pd.merge(data_subset, glo_subset, on=columns_to_compare, how='inner')
>>> compare_sametrans = pd.DataFrame({'chrom':sametrans['chrom'],'chromStart':sametrans['chromStart'],'chromEnd':sametrans['chromEnd'],'strand':sametrans['strand'],'GLORI_score':glo['score'],'m6A_frequency':data['frequency']})
>>> print(compare_sametrans.head())
  chrom  chromStart  chromEnd strand  GLORI_score  m6A_frequency
0     1    185037.0  185038.0      -         31.1            0.0
1     1    185232.0  185233.0      -         15.4            0.0
2     1    826938.0  826939.0      -         25.7            0.0
3     1    841506.0  841507.0      +         25.7            0.0
4     1    841558.0  841559.0      +         18.3            0.0
>>> correlation = compare_sametrans[['GLORI_score', 'm6A_frequency']].corr()
>>> print(correlation)
               GLORI_score  m6A_frequency
GLORI_score       1.000000       0.003066
m6A_frequency     0.003066       1.000000
>>> sns.set(style="white")
>>> plt.figure(figsize=(8, 6))
<Figure size 800x600 with 0 Axes>
>>> sns.heatmap(correlation, annot=True, cmap="coolwarm", fmt=".2f", linewidths=0.5)
<Axes: >
>>> plt.title("Correlation between GLORI_score and m6A_frequency")
Text(0.5, 1.0, 'Correlation between GLORI_score and m6A_frequency')
>>> plt.savefig('heatmap.png')


