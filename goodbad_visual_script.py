>>> import pandas as pd
>>> import matplotlib.pyplot as plt
>>> import seaborn as sns
>>> df = pd.read_csv("grouped_sites_updated.csv", sep='\t',low_memory=False)
>>> print(df.head(10))
  chrom  chromStart  chromEnd strand  GLORI_score  m6A_frequency  relative_difference group
0     1      185037    185038      -         69.3          61.11                0.118  good
1     1      185232    185233      -         34.4          12.09                0.649   bad
2     1      826938    826939      -         65.8          57.14                0.132  good
3     1      841506    841507      +         73.7          92.31                0.202   bad
4     1      841558    841559      +         74.5          91.67                0.187  good
5     1      841599    841600      +         80.6          92.31                0.127  good
6     1      841882    841883      +         30.6          42.86                0.286   bad
7     1      842526    842527      +         95.8         100.00                0.042  good
8     1      842532    842533      +         55.4          50.00                0.097  good
9     1      842559    842560      +         26.1          81.82                0.681   bad
>>> print(df.columns.tolist())
['chrom', 'chromStart', 'chromEnd', 'strand', 'GLORI_score', 'm6A_frequency', 'relative_difference', 'group']
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
       34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47]), [Text(0, 0, '(1, +)'), Text(1, 0, '(1, -)'), Text(2, 0, '(10, +)'), Text(3, 0, '(10, -)'), Text(4, 0, '(11, +)'), Text(5, 0, '(11, -)'), Text(6, 0, '(12, +)'), Text(7, 0, '(12, -)'), Text(8, 0, '(13, +)'), Text(9, 0, '(13, -)'), Text(10, 0, '(14, +)'), Text(11, 0, '(14, -)'), Text(12, 0, '(15, +)'), Text(13, 0, '(15, -)'), Text(14, 0, '(16, +)'), Text(15, 0, '(16, -)'), Text(16, 0, '(17, +)'), Text(17, 0, '(17, -)'), Text(18, 0, '(18, +)'), Text(19, 0, '(18, -)'), Text(20, 0, '(19, +)'), Text(21, 0, '(19, -)'), Text(22, 0, '(2, +)'), Text(23, 0, '(2, -)'), Text(24, 0, '(20, +)'), Text(25, 0, '(20, -)'), Text(26, 0, '(21, +)'), Text(27, 0, '(21, -)'), Text(28, 0, '(22, +)'), Text(29, 0, '(22, -)'), Text(30, 0, '(3, +)'), Text(31, 0, '(3, -)'), Text(32, 0, '(4, +)'), Text(33, 0, '(4, -)'), Text(34, 0, '(5, +)'), Text(35, 0, '(5, -)'), Text(36, 0, '(6, +)'), Text(37, 0, '(6, -)'), Text(38, 0, '(7, +)'), Text(39, 0, '(7, -)'), Text(40, 0, '(8, +)'), Text(41, 0, '(8, -)'), Text(42, 0, '(9, +)'), Text(43, 0, '(9, -)'), Text(44, 0, '(MT, +)'), Text(45, 0, '(MT, -)'), Text(46, 0, '(X, +)'), Text(47, 0, '(X, -)')])
>>> plt.legend(title="Group and Strand")
<matplotlib.legend.Legend object at 0x7f303bda59d0>
>>> plt.tight_layout()
>>> plt.savefig("chromosome_strand_distribution.png")
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
<matplotlib.legend.Legend object at 0x7f3038a83a10>
>>> plt.tight_layout()
<stdin>:1: UserWarning: Creating legend with loc="best" can be slow with large amounts of data.
>>> plt.savefig("chromosome_strand_scatter.png")
>>> # pie chart
>>> strand_counts = df.groupby(['strand', 'group']).size().unstack(fill_value=0)
>>> strand_counts.plot(kind='pie', subplots=True, autopct='%1.1f%%', figsize=(12, 6))
array([<Axes: ylabel='bad'>, <Axes: ylabel='good'>], dtype=object)
>>> plt.title("Strand Distribution of Good and Bad Predictions")
Text(0.5, 1.0, 'Strand Distribution of Good and Bad Predictions')
>>> plt.tight_layout()
>>> plt.savefig("strand_distribution_pie.png")
>>> exit()
