#seqanno

Retrieves read counts originating at a given genomic sequence and allows further characterization. 

##Usage

###bam2bedgraph

Convert bam to bedgraph format. Options for unstranded or stranded counts.

```python
python seqanno.py bedgraph bps.bam > bps.bedgraph
```

###search

Search for a genomic sequence, report read counts and relative location within a
gene.

```python
python seqanno.py search bps.bedgraph bps.fa bps.gff TA > annotated_sites.bed
```

###genestat

Get gene level count statistics.

```python
python seqanno.py genestat annotated_sites.bed 5 5,80 > gene_level_counts.txt
```

###compare

Retrieve genes unique to file _A_.

```python
python seqanno.py compare a.test b.test aonly 5 > a_unique.txt
```

###uniprot

Append Uniprot data onto tab delimited text file. Specify the column to match to
the Uniprot 'Gene names' column.

```python
python seqanno.py uniprot gene_stats.txt bps_uniprot.txt 1 > uniprot_columns_appended.txt

python seqanno.py uniprot annotated_sites.bed bps_uniprot.txt 5 > annotated_sites_with_uniprot.txt
```