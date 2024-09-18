# Useful scripts for bioinformatic tasks
These are some scripts that I use multiple times. They are not perfect or efficient but they work. 

__1) check_rnaseq_data.py__
This script uses RSeQC (https://doi.org/10.1093/bioinformatics/bts356) to determine whether a set of RNAseq data is stranded.

__2) check_te_overlap_sequences.py__
Check the TE overlap between a TE annotation (bedfile) and a sequence file (bedfile). The output will be two files: one containing the relative overlap with TE and the other containing information about the TE types.

__3) gff2bed.py__
Transform a gff file to a bedfile. 
