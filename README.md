# deswomanUTIL
Some scripts to deal with the DESwoMAN output.

> [!IMPORTANT]
> This script has not been extensively tested so use at your own risk and double check the results. 

## 1) filterDESwoMAN.py
Filter the output of [DESwoMAN](https://github.com/AnnaGrBio/DESWOMAN) to get a dataset of high confidence _de novo_ originated neORFs.
It is important that you did the following steps before running:
1) You ran DESwoMAN on one or multiple species. The output for each species needs to be in a separate folder with the full name (e.g. DmelZI or DsubFAL). This is also true when you only ran it for one species!
2) You have a newick tree with all species (in- and outgroups) in your analysis (even if you also have populations this tree needs to contain 4 letter species ids). It is important that all internal nodes are named. Also no polytonies are allowed as this will lead to errors!
3) Your samples have a four letter species id optionally followed by a sample id (if applicable) e.g. Dmel (for _Drosophila melanogaster_) or Hsap (for _Homo sapiens_).

- `--te_check` Blast against an TE database  (default = false)
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
- `--accepted_mutations` Comma separated list of all mutations you want to count as noncoding (all possibilities are: start, stop, prematurestop and transcription (complete, partial and/or antisense))
- `--frameshift_score` Maximum percent of frame conservation above which it will be considered coding

__Usage:__
```python3 filterDESwoMAN.py --te_db Database.fa --ortho Orthofinder/OrthoFinder/Results/Orthogroups.txt --deswoman ./ --out FilteredResults```


## 2) getGFFfromOut.py
Generate a gff file with all neORF candidates from the DESwoMAN info file. 

- `--deswoman` Path to the deswoman information file
- `--gtf` Path to the corresponding transcriptome assembly (gtf)
- `--outname` Name of the output file (default: DESwoMAN)

__Usage:__
```python3 getGFFfromOut.py [-h] --deswoman PATH/TO/DESWOMAN/info_file.txt --gtf Path/To/GTF --outname Path/To/OUTNAME```