>>> import pandas as pd
>>> df = pd.read_csv("BIDpsi_merged_data.bed", sep='\t', low_memory=False)
>>> print(df.head(10))
  chrom  chromStart  chromEnd strand  BID_score  psi_frequency
0     1     1045780   1045781      +        7.7           3.03
1     1     1051519   1051520      +       23.4          61.54
2     1     1082417   1082418      -       63.9          20.00
3     1     1254956   1254957      -       13.1          14.09
4     1     1336109   1336110      -        9.3           1.00
5     1     2589932   2589933      +       10.3           3.17
6     1     6185969   6185970      -       18.4          46.09
7     1     6186767   6186768      -       10.2          22.72
8     1     6222864   6222865      -       26.7          30.00
9     1     6623606   6623607      +       51.2          51.95
>>> print(df.columns.tolist())
['chrom', 'chromStart', 'chromEnd', 'strand', 'BID_score', 'psi_frequency']
>>> df['difference'] = abs(df['BID_score'] - df['psi_frequency'])
>>> df['relative_difference'] = df['difference'] / df[['BID_score', 'psi_frequency']].max(axis=1)
>>> df['relative_difference'] = df['relative_difference'].round(3)
>>> threshold = 0.2
>>> df["group"] = df["relative_difference"].apply(lambda x: "good" if x <= threshold else "bad")
>>> df.to_csv("BIDpsi_grouped_sites.csv", index=False)
>>> group = pd.read_csv("BIDpsi_grouped_sites.csv", sep='\t', low_memory=False)
>>> print(df.columns.tolist())
['chrom', 'chromStart', 'chromEnd', 'strand', 'BID_score', 'psi_frequency', 'difference', 'relative_difference', 'group']
>>> df = df.drop(columns=['difference'])
>>> df.to_csv("BIDpsigrouped_sites_updated.csv", index=False, sep='\t')
>>> exit()
