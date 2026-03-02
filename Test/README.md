# Test dataset

If you are unsure how to structure your data to run the scripts you can have a look here. For each strategy (1 and 2), one exemplary output folder is provided.

-------------------
The following genomes were used to create the DESwoMAN dataset:<br>

**Strategy1:**<br>
- Queries:  _Drosophila melanogaster_: FI and ZI genomes and transcriptomes ([downloaded from Grandchamp et al. (2023)](https://zenodo.org/records/7322757))  
- Targets: 
    - _Drosophila simulans_ (GCF_016746395.2)
    - _Drosophila suzukii_ (GCF_037355615.1)
    - _Drosophila ananassae_ (GCF_017639315.1)

**Strategy2**<br>
- Query: _Drosophila melanogaster_ reference genome (GCF_000001215.4)
- Transcriptomes: RNAseq from FI and ZI was downloaded from NCBI (Bioproject: PRJNA929424, more detail on data in [Grandchamp et al. (2023)](https://zenodo.org/records/7322757)) and mapped to the reference genome. Then the transcriptomes (for each line) were extracted with stringtie.

DESwoMAN was run for both strategies using the following parameters:
- Intergenic ORFs only
- TPM >= 0.5 (+ Exclude undirectional transcripts ".")
- Protein search (--more-sensitive option) against _Drosophila_ and Diptera from Ensembl Metazoa (database not uploaded due to size, available upon request)
- Synteny window = 4 (no reciprocal blast) - _Strategy 1 only_

Additional data downloaded:
- [Flybase TE database](https://flybase.org/downloads/bulkdata) (accessed July, 2025)
- ncRNA database: [Drosophila simulans ncRNA](https://metazoa.ensembl.org/Drosophila_simulans_gca016746395v2rs/Info/Index) (Assembly: GCF_016746395.2)

There are several files and folders in the test dataset.
- `DESwoMAN_Strat1:` Contains the DESwoMAN Strategy 1 output for both query lines (FI and ZI) (Intermediate output removed due to size)
- `DESwoMAN_Strat2:` Contains the DESwoMAN Strategy 2 output for both query lines (FI and ZI) (Intermediate output removed due to size)
- `Transcriptomes_Strat1:` Contains the query transcriptome gtf files to run Strategy 1(downloaded from zenodo)
- `Transcriptomes_Strat2:` Contains the query transcriptome gtf files to run Strategy (mapped to reference genome + assembled with stringtie)
- `FlyBase_TE.fa:` Flybase TE database
- `NcRNA:` Noncoding RNA database (= for testing purposes only ncRNA of _Drosophila simulans_)
- `Orthogroups.txt:` Orthogroups generated with Orthofinder using the script [runOrthofinder.py](https://github.com/MarieLebh/deswomanUTIL/blob/main/runOrthfinder.py)
- `Query_names.txt:` File with all query species (= species file for all provided scripts)
- `Tree.nwk:` A newick tree of all query and target species. Important: All internal nodes need an Id and no polytonies are allowed. 