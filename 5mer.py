>>> from Bio import SeqIO
>>> genome_reference = "/biodb/genomes/homo_sapiens/GRCh38_102/GRCh38_102.fa"
>>> genome_dict = {}
>>> for record in SeqIO.parse(genome_reference, "fasta"):
...     genome_dict[record.id] = record.seq
... 
>>> from Bio.Seq import Seq
>>> def reverse_complement(sequence):
...     return str(Seq(sequence).reverse_complement())
... 
>>> import pandas as pd
>>> df = pd.read_csv("grouped_sites_updated.csv", sep='\t', low_memory=False)
>>> print(df.columns.tolist())
['chrom', 'chromStart', 'chromEnd', 'strand', 'GLORI_score', 'm6A_frequency', 'relative_difference', 'group']
>>> def get_5mer_context(row):
...     chrom = row['chrom']
...     start = row['chromStart']
...     strand = row['strand']
...     if chrom not in genome_dict:
...             return "NA"
...     chrom_sequence = genome_dict[chrom]
...     if start < 2 or start >= len(chrom_sequence) - 2:
...             return "NA"
...     context = chrom_sequence[start - 2:start + 3]
...     if strand == '-':
...             context = reverse_complement(context)
...     return context
... 
>>> df['5mer_context'] = df.apply(get_5mer_context, axis=1)
>>> df.to_csv("grouped_sites_with_5mer.csv", sep='\t', index=False)
>>> fivemer = pd.read_csv("grouped_sites_with_5mer.csv", sep='\t', low_memory=False)
>>> print(fivemer.head(10))
  chrom  chromStart  chromEnd strand  GLORI_score  m6A_frequency  relative_difference group 5mer_context
0     1      185037    185038      -         69.3          61.11                0.118  good        GGACT
1     1      185232    185233      -         34.4          12.09                0.649   bad        GGACG
2     1      826938    826939      -         65.8          57.14                0.132  good        TGACT
3     1      841506    841507      +         73.7          92.31                0.202   bad        AGACT
4     1      841558    841559      +         74.5          91.67                0.187  good        GAACT
5     1      841599    841600      +         80.6          92.31                0.127  good        GGACT
6     1      841882    841883      +         30.6          42.86                0.286   bad        CAACC
7     1      842526    842527      +         95.8         100.00                0.042  good        GGACT
8     1      842532    842533      +         55.4          50.00                0.097  good        TGACT
9     1      842559    842560      +         26.1          81.82                0.681   bad        AGACA
>>> exit()
