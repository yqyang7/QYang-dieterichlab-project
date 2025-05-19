>>> import pandas as pd
>>> df = pd.read_csv("file_merged_data.bed", sep='\t', low_memory=False)
>>> print(df.head(10))
  chrom  chromStart  chromEnd strand  GLORI_score  m6A_frequency
0     1      185037    185038      -         69.3          61.11
1     1      185232    185233      -         34.4          12.09
2     1      826938    826939      -         65.8          57.14
3     1      841506    841507      +         73.7          92.31
4     1      841558    841559      +         74.5          91.67
5     1      841599    841600      +         80.6          92.31
6     1      841882    841883      +         30.6          42.86
7     1      842526    842527      +         95.8         100.00
8     1      842532    842533      +         55.4          50.00
9     1      842559    842560      +         26.1          81.82
>>> df['difference'] = abs(df['GLORI_score'] - df['m6A_frequency'])
>>> df['relative_difference'] = df['difference'] / df[['GLORI_score', 'm6A_frequency']].max(axis=1)
>>> df['relative_difference'] = df['relative_difference'].round(3)
>>> threshold = 0.2
>>> df["group"] = df["relative_difference"].apply(lambda x: "good" if x <= threshold else "bad")
>>> df.to_csv("grouped_sites.csv", index=False)
>>> group = pd.read_csv("grouped_sites.csv", sep='\t', low_memory=False)
>>> df = pd.read_csv("grouped_sites.csv", sep=',', low_memory=False)
>>> print(df.columns.tolist())
['chrom', 'chromStart', 'chromEnd', 'strand', 'GLORI_score', 'm6A_frequency', 'difference', 'relative_difference', 'group']
>>> df = df.drop(columns=['difference'])
>>> df.to_csv("grouped_sites_updated.csv", index=False, sep='\t')
>>> exit()
