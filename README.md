# Useful scripts for bioinformatic tasks
Here are some scripts I wrote for tasks I needed to run multiple times.

> [!IMPORTANT]
> This script has not been extensively tested so use at your own risk and double check the results. 

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

## 4) filterDESwoMAN.py
Filter the output of [DESwoMAN](https://github.com/AnnaGrBio/DESWOMAN) to get a dataset of high confidence _de novo_ originated neORFs.
It is important that you did the following steps before running:
1) You ran DESwoMAN on one or multiple species.
2) You ran Orthofinder for all query species in your DESwoMAN runs and used the resulting neORFs to get Orthogroups
3) You have a newick tree with all species (in- and outgroups) in your analysis
4) Your samples have a four letter species id optionally followed by a sample id (if applicable) e.g. Dmel (for _Drosophila melanogaster_) or Hsap (for _Homo sapiens_). 

- `--te_db` Path to a fasta TE database of your organism (str)
- `--te_cov` Blast coverage for TE search (float, default = 80)
- `--te_idt` Percent identity for TE search (float, default = 80)
- `--te_eval` Evalue for TE search (float, default = 0.001)
- `--ortho` Path to an Orthogroup file ("Orthogroups.txt", str) generated with Orthofinder and including the filename
- `--deswoman` Path to the DESwoMAN folder. Its necessary that you run DESwoMAN for all species.
- `--rna_check` Blast against an rna database (e.g. cdna/ncrna) (default = false)
- `--tr_db`  Path to a transcript database (e.g. cdna/ncrna, str)
- `--tr_cov` Blast coverage for transcript search (float, default = 80)
- `--tr_idt` Percent identity for transcript search (float, default = 80)
- `--tr_eval` Evalue for transcript search (float, default = 0.001)
- `--tr_strand` Strand for the transcript search (string, default = "plus")
- `--out` Name of the output folder to store the filtered data
- `--tree` Newick tree of all samples analyzed (needs internal node ids, str)
- `--species_file` Text file with all samples listed (one per line, str)

__Usage:__
```python3 GetHighConfidenceNeORFs.py --te_db Database.fa --ortho Orthofinder/OrthoFinder/Results/Orthogroups.txt --deswoman ./ --out FilteredResults```
