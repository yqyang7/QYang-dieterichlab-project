>>> import pandas as pd
>>> import matplotlib
>>> matplotlib.use('Agg')
>>> import seaborn as sns
>>> import matplotlib.pyplot as plt
>>> psi = pd.read_csv("/prj/Qingyang/psi_cov10.bed", sep=r'\s+', low_memory=False)
>>> print(psi.head())
  chrom  chromStart  chromEnd   name  score strand  thickStart  thickEnd  itemRgb  coverage  frequency
0     1       14384     14385  17802     11      -       14384     14385  255,0,0        11       9.09
1     1       14387     14388  17802     15      -       14387     14388  255,0,0        15       0.00
2     1       14390     14391  17802     10      -       14390     14391  255,0,0        10      30.00
3     1       14393     14394  17802     14      -       14393     14394  255,0,0        14       7.14
4     1       14395     14396  17802     12      -       14395     14396  255,0,0        12       0.00
>>> BID = pd.read_csv("/prj/Qingyang/BID_seq_HEK293T.bed", sep=r'\s+', low_memory=False)
>>> print(BID.head())
  chrom  chromStart  chromEnd name  score strand ref5mer  BID-Seq
0     1     1045780   1045781  psi    7.7      +   GATGA     12.4
1     1     1051519   1051520  psi   23.4      +   CTTTG     73.3
2     1     1082417   1082418  psi   63.9      -   GTTCC     90.5
3     1     1254956   1254957  psi   13.1      -   TGTAG     45.5
4     1     1336109   1336110  psi    9.3      -   GGTGG     12.4
>>> columns_to_compare = ["chrom","chromStart","chromEnd","strand"]
>>> psi_subset = psi[columns_to_compare]
>>> BID_subset = BID[columns_to_compare]
>>> sametrans = pd.merge(psi_subset, BID_subset, on=columns_to_compare, how='inner')
>>> compare_sametrans = pd.DataFrame({'chrom':sametrans['chrom'],'chromStart':sametrans['chromStart'],'chromEnd':sametrans['chromEnd'],'strand':sametrans['strand'],'psi_frequency':psi['frequency'],'BID_score':BID['score']})
>>> print(compare_sametrans.head())
  chrom  chromStart   chromEnd strand  psi_frequency  BID_score
0     1   1045780.0  1045781.0      +           9.09        7.7
1     1   1051519.0  1051520.0      +           0.00       23.4
2     1   1082417.0  1082418.0      -          30.00       63.9
3     1   1254956.0  1254957.0      -           7.14       13.1
4     1   1336109.0  1336110.0      -           0.00        9.3
>>> correlation = compare_sametrans[['BID_score', 'psi_frequency']].corr()
>>> print(correlation)
               BID_score  psi_frequency
BID_score       1.000000       0.063347
psi_frequency   0.063347       1.000000
>>> plt.figure(figsize=(8, 6))
<Figure size 800x600 with 0 Axes>
>>> plt.scatter(BID_score, psi_frequency, alpha=0.6)
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
NameError: name 'BID_score' is not defined
>>> plt.scatter(compare_sametrans['BID_score'],compare_sametrans['psi_frequency'], alpha = 0.6)
<matplotlib.collections.PathCollection object at 0x7f1783753110>
>>> plt.title('Scatterplot of BID_score vs psi_frequency')
Text(0.5, 1.0, 'Scatterplot of BID_score vs psi_frequency')
>>> plt.xlabel('BID_score')
Text(0.5, 0, 'BID_score')
>>> plt.ylabel('psi_frequency')
Text(0, 0.5, 'psi_frequency')
>>> plt.grid(True)
>>> plt.show()

