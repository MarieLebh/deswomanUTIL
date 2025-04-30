# Useful scripts for bioinformatic tasks
Here are some scripts I wrote for tasks I needed to run multiple times.

## 1) check_rnaseq_data.py
This script uses [RSeQC](https://doi.org/10.1093/bioinformatics/bts356) to determine whether a set of RNAseq data is stranded.

## 2) Fai_converter.py
Convert a genome index file (.fai) to a chromosome length bed3file (ChromName 0 Chrom Length).

- `--inp` Path to the index file including the filename (e.g. /home/user/Genome.fa.fai)
- `--out` Path to the output file including the filename (e.g. /home/user/Genome.bed) - _optional_

__Usage:__
```python3 Fai_converter.py -in PathToGenome.fai -out PathToOutput.bed (default = PathToGenome_ChromLengths.bed)```

## 3) gff2bed.py
Transform a gff/gtf file to a bedfile (bed6). 

- `--inp` Path to the gtf/gff file including the filename (e.g. /home/user/genes.gff)
- `--out` Path to the output file excluding the filename (e.g. /home/user/genes)
- `--type` Type (e.g. gene, transcript, mRNA, ...) to be included in the bedfile (only 1 possible)  - _optional, default: mRNA_

 __Usage:__
```python3 gff2bed.py --inp GFFFile (--out Outname (no suffix)) (--type Type (e.g. gene/mRNA/transcript, ...))```
