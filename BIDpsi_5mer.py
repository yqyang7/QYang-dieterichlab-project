>>> from Bio import SeqIO
>>> genome_reference = "/biodb/genomes/homo_sapiens/GRCh38_102/GRCh38_102.fa"
>>> genome_dict = {}
>>> for record in SeqIO.parse(genome_reference, "fasta"):
...     genome_dict[record.id] = record.seq
... 

>>> 
>>> from Bio.Seq import Seq
>>> def reverse_complement(sequence):
...     return str(Seq(sequence).reverse_complement())
... 
>>> import pandas as pd
>>> df = pd.read_csv("BIDpsigrouped_sites_updated.csv", sep='\t', low_memory=False)
>>> print(df.columns.tolist())
['chrom', 'chromStart', 'chromEnd', 'strand', 'BID_score', 'psi_frequency', 'relative_difference', 'group']
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
>>> df.to_csv("BIDpsi_5mer.csv", sep='\t', index=False)
>>> fivemer = pd.read_csv("BIDpsi_5mer.csv", sep='\t', low_memory=False)
>>> print(fivemer.head(10))
  chrom  chromStart  chromEnd strand  BID_score  psi_frequency  relative_difference group 5mer_context
0     1     1045780   1045781      +        7.7           3.03                0.606   bad        GATGA
1     1     1051519   1051520      +       23.4          61.54                0.620   bad        CTTTG
2     1     1082417   1082418      -       63.9          20.00                0.687   bad        GTTCC
3     1     1254956   1254957      -       13.1          14.09                0.070  good        TGTAG
4     1     1336109   1336110      -        9.3           1.00                0.892   bad        GGTGG
5     1     2589932   2589933      +       10.3           3.17                0.692   bad        GGTGG
6     1     6185969   6185970      -       18.4          46.09                0.601   bad        TGTAG
7     1     6186767   6186768      -       10.2          22.72                0.551   bad        CGTGA
8     1     6222864   6222865      -       26.7          30.00                0.110  good        GTTAG
9     1     6623606   6623607      +       51.2          51.95                0.014  good        GTTCG
>>> exit()
