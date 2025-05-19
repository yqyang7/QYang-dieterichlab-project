>>> import pandas as pd
>>> import matplotlib
>>> matplotlib.use('Agg')
>>> import seaborn as sns
>>> import matplotlib.pyplot as plt
>>> psi = pd.read_csv("/prj/Qingyang/psi_cov10.bed", sep=r'\s+', low_memory=False)
>>> print(psi.head(10))
  chrom  chromStart  chromEnd   name  score strand  thickStart  thickEnd  itemRgb  coverage  frequency
0     1       14384     14385  17802     11      -       14384     14385  255,0,0        11       9.09
1     1       14387     14388  17802     15      -       14387     14388  255,0,0        15       0.00
2     1       14390     14391  17802     10      -       14390     14391  255,0,0        10      30.00
3     1       14393     14394  17802     14      -       14393     14394  255,0,0        14       7.14
4     1       14395     14396  17802     12      -       14395     14396  255,0,0        12       0.00
5     1       14412     14413  17802    153      -       14412     14413  255,0,0       153       1.31
6     1       14420     14421  17802    130      -       14420     14421  255,0,0       130       8.46
7     1       14424     14425  17802    152      -       14424     14425  255,0,0       152       1.32
8     1       14435     14436  17802     16      -       14435     14436  255,0,0        16       0.00
9     1       14446     14447  17802    132      -       14446     14447  255,0,0       132       9.09
>>> BID = pd.read_csv("/prj/Qingyang/BID_seq_HEK293T.bed", sep=r'\s+', low_memory=False)
>>> print(BID.head(10))
  chrom  chromStart  chromEnd name  score strand ref5mer  BID-Seq
0     1     1045780   1045781  psi    7.7      +   GATGA     12.4
1     1     1051519   1051520  psi   23.4      +   CTTTG     73.3
2     1     1082417   1082418  psi   63.9      -   GTTCC     90.5
3     1     1254956   1254957  psi   13.1      -   TGTAG     45.5
4     1     1336109   1336110  psi    9.3      -   GGTGG     12.4
5     1     2589932   2589933  psi   10.3      +   GGTGG     15.4
6     1     6185969   6185970  psi   18.4      -   TGTAG     59.0
7     1     6186767   6186768  psi   10.2      -   CGTGA     21.3
8     1     6222864   6222865  psi   26.7      -   GTTAG     53.5
9     1     6623606   6623607  psi   51.2      +   GTTCG     81.8
>>> merged_data = pd.merge(psi, BID, on=['chrom', 'chromStart', 'chromEnd', 'strand'],how='inner',suffixes=('_psi', '_BID'))
>>> result = merged_data[['chrom', 'chromStart', 'chromEnd', 'strand', 'score_BID', 'frequency']]
>>> result.columns = ['chrom', 'chromStart', 'chromEnd', 'strand', 'BID_score', 'psi_frequency']
>>> result.to_csv("BIDpsi_merged_data.bed", index=False, sep='\t')
>>> print(result)
    chrom  chromStart   chromEnd strand  BID_score  psi_frequency
0       1     1045780    1045781      +        7.7           3.03
1       1     1051519    1051520      +       23.4          61.54
2       1     1082417    1082418      -       63.9          20.00
3       1     1254956    1254957      -       13.1          14.09
4       1     1336109    1336110      -        9.3           1.00
..    ...         ...        ...    ...        ...            ...
494     X   135022169  135022170      -       21.2           0.00
495     X   153700679  153700680      -       31.4           2.36
496     X   153703028  153703029      -       20.5          23.45
497     X   154353116  154353117      -        8.9           0.58
498     X   154765960  154765961      +       52.6          79.78

[499 rows x 6 columns]
>>> correlation = result['BID_score'].corr(result['psi_frequency'])
>>> print(f"\correlation between BID_score and psi_frequency：{correlation:.4f}")
\correlation between BID_score and psi_frequency：0.8278
>>> plt.figure(figsize=(8, 6))
<Figure size 800x600 with 0 Axes>
>>> sns.scatterplot(x='BID_score',y='psi_frequency',data=result,color="blue",alpha=0.7,edgecolor=None)
<Axes: xlabel='BID_score', ylabel='psi_frequency'>
>>> plt.title("Scatterplot of BID_score vs psi_frequency", fontsize=14)
Text(0.5, 1.0, 'Scatterplot of BID_score vs psi_frequency')
>>> plt.xlabel("BID_score", fontsize=12)
Text(0.5, 0, 'BID_score')
>>> plt.ylabel("psi_frequency", fontsize=12)
Text(0, 0.5, 'psi_frequency')
>>> sns.regplot(x='BID_score',y='psi_frequency',data=result,scatter=False,color="red",line_kws={'linewidth': 2})
<Axes: title={'center': 'Scatterplot of BID_score vs psi_frequency'}, xlabel='BID_score', ylabel='psi_frequency'>
>>> plt.tight_layout()
>>> plt.show()
>>> plt.savefig('scatterplot of BID_score vs psi_frequency.png')
>>> b = result['BID_score']
>>> p = result['psi_frequency']
>>> plt.figure(figsize=(10, 8))
<Figure size 1000x800 with 0 Axes>
>>> plt.hexbin(b, p, gridsize=50, cmap='coolwarm') 
<matplotlib.collections.PolyCollection object at 0x7f84adb4d1d0>
>>> plt.colorbar(label='Density')
<matplotlib.colorbar.Colorbar object at 0x7f84adbd6410>
>>> plt.xlabel('BID_score')
Text(0.5, 0, 'BID_score')
>>> plt.ylabel('psi_frequency')
Text(0, 0.5, 'psi_frequency')
>>> plt.title('Hexbin Heatmap: BID_score vs psi_frequency')
Text(0.5, 1.0, 'Hexbin Heatmap: BID_score vs psi_frequency')
>>> plt.show()
>>> plt.savefig('BIDpsi heatmap_new.png')
>>> exit()
