>>> import pandas as pd
>>> import matplotlib.pyplot as plt
>>> import seaborn as sns
>>> df = pd.read_csv("BIDpsigrouped_sites_updated.csv", sep='\t',low_memory=False)
>>> print(df.head(10))
  chrom  chromStart  chromEnd strand  BID_score  psi_frequency  relative_difference group
0     1     1045780   1045781      +        7.7           3.03                0.606   bad
1     1     1051519   1051520      +       23.4          61.54                0.620   bad
2     1     1082417   1082418      -       63.9          20.00                0.687   bad
3     1     1254956   1254957      -       13.1          14.09                0.070  good
4     1     1336109   1336110      -        9.3           1.00                0.892   bad
5     1     2589932   2589933      +       10.3           3.17                0.692   bad
6     1     6185969   6185970      -       18.4          46.09                0.601   bad
7     1     6186767   6186768      -       10.2          22.72                0.551   bad
8     1     6222864   6222865      -       26.7          30.00                0.110  good
9     1     6623606   6623607      +       51.2          51.95                0.014  good
>>> print(df.columns.tolist())
['chrom', 'chromStart', 'chromEnd', 'strand', 'BID_score', 'psi_frequency', 'relative_difference', 'group']
>>> chrom_counts = df.groupby(['chrom', 'group']).size().unstack(fill_value=0)
>>> # Histogram (chromosome and strand statistics)
>>> grouped_data = df.groupby(['chrom', 'strand', 'group']).size().unstack(fill_value=0)
>>> grouped_data.plot(kind='bar', stacked=True, figsize=(12, 6), colormap='coolwarm')
<Axes: xlabel='chrom,strand'>
>>> plt.title("Good and Bad Predictions Across Chromosomes and Strands")
Text(0.5, 1.0, 'Good and Bad Predictions Across Chromosomes and Strands')
>>> plt.ylabel("Count")
Text(0, 0.5, 'Count')
>>> plt.xlabel("Chromosome")
Text(0.5, 0, 'Chromosome')
>>> plt.xticks(rotation=45)
(array([ 0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16,
       17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33,
       34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45]), [Text(0, 0, '(1, +)'), Text(1, 0, '(1, -)'), Text(2, 0, '(10, +)'), Text(3, 0, '(10, -)'), Text(4, 0, '(11, +)'), Text(5, 0, '(11, -)'), Text(6, 0, '(12, +)'), Text(7, 0, '(12, -)'), Text(8, 0, '(13, +)'), Text(9, 0, '(13, -)'), Text(10, 0, '(14, +)'), Text(11, 0, '(14, -)'), Text(12, 0, '(15, +)'), Text(13, 0, '(15, -)'), Text(14, 0, '(16, +)'), Text(15, 0, '(16, -)'), Text(16, 0, '(17, +)'), Text(17, 0, '(17, -)'), Text(18, 0, '(18, +)'), Text(19, 0, '(18, -)'), Text(20, 0, '(19, +)'), Text(21, 0, '(19, -)'), Text(22, 0, '(2, +)'), Text(23, 0, '(2, -)'), Text(24, 0, '(20, +)'), Text(25, 0, '(20, -)'), Text(26, 0, '(21, +)'), Text(27, 0, '(21, -)'), Text(28, 0, '(22, +)'), Text(29, 0, '(22, -)'), Text(30, 0, '(3, +)'), Text(31, 0, '(3, -)'), Text(32, 0, '(4, +)'), Text(33, 0, '(4, -)'), Text(34, 0, '(5, +)'), Text(35, 0, '(5, -)'), Text(36, 0, '(6, +)'), Text(37, 0, '(6, -)'), Text(38, 0, '(7, +)'), Text(39, 0, '(7, -)'), Text(40, 0, '(8, +)'), Text(41, 0, '(8, -)'), Text(42, 0, '(9, +)'), Text(43, 0, '(9, -)'), Text(44, 0, '(X, +)'), Text(45, 0, '(X, -)')])
>>> plt.legend(title="Group and Strand")
<matplotlib.legend.Legend object at 0x7f5e7f8378d0>
>>> plt.tight_layout()
>>> plt.savefig("BIDpsi_chrom_str_distribution.png")
>>> # Scatterplot (chromosome region and strand direction)
>>> plt.figure(figsize=(12, 6))
<Figure size 1200x600 with 0 Axes>
>>> sns.scatterplot(data=df,x='chromStart',y='chrom',hue='group',style='strand',palette='coolwarm',s=50)
<Axes: xlabel='chromStart', ylabel='chrom'>
>>> plt.title("Good and Bad Predictions Across Chromosome Regions and Strands")
Text(0.5, 1.0, 'Good and Bad Predictions Across Chromosome Regions and Strands')
>>> plt.xlabel("Chromosome Start Position")
Text(0.5, 0, 'Chromosome Start Position')
>>> plt.ylabel("Chromosome")
Text(0, 0.5, 'Chromosome')
>>> plt.legend(title="Group and Strand")
<matplotlib.legend.Legend object at 0x7f5e7d368fd0>
>>> plt.tight_layout()
>>> plt.savefig("BIDpsi_chrom_str_scatter.png")
>>> # pie chart
>>> strand_counts = df.groupby(['strand', 'group']).size().unstack(fill_value=0)
>>> strand_counts.plot(kind='pie', subplots=True, autopct='%1.1f%%', figsize=(12, 6))
array([<Axes: ylabel='bad'>, <Axes: ylabel='good'>], dtype=object)
>>> plt.title("Strand Distribution of Good and Bad Predictions")
Text(0.5, 1.0, 'Strand Distribution of Good and Bad Predictions')
>>> plt.tight_layout()
>>> plt.savefig("BIDpsi-str_distri_pie.png")
>>> exit()
