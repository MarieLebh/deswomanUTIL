# Useful scripts for bioinformatic tasks
Here are some scripts I wrote for tasks I needed to run multiple times.

__1) check_rnaseq_data.py__
This script uses [RSeQC](https://doi.org/10.1093/bioinformatics/bts356) to determine whether a set of RNAseq data is stranded.

__2) Fai_converter.py__
Convert a genome index file (.fai) to a chromosome length bed3file (ChromName 0 Chrom Length).

```python3 Fai_converter.py -in PathToGenome.fai -out PathToOutput.bed (default = PathToGenome_ChromLengths.bed)```

__3) gff2bed.py__
Transform a gff file to a bedfile. 

```python3 gff2bed.py --inp GFFFile (--out Outname (no suffix)) (--type Type (e.g. gene/mRNA/transcript, ...))```
